version 1.0
    
import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

# cohort set
workflow inferPlatform {
    input {
        Array[File] vcf_files
        Array[String] cohort_prefixes
        File interval_list

        String cohort_set_id
        String bucket_id
        String hail_docker
        String genome_build='GRCh38'

        Float binarization_threshold=0.25
        Int n_pcs=10
        Int hdbscan_min_cluster_size=0
        Int hdbscan_min_samples=50
        RuntimeAttr? runtime_attr_get_call_rate_mt
        RuntimeAttr? runtime_attr_infer_platform
    }

    scatter (pair in zip(vcf_files, cohort_prefixes)) {
        File vcf_file = pair.left
        String cohort_prefix = pair.right
        call getCallRateMT {
            input:
            vcf_file=vcf_file,
            interval_list=interval_list,
            cohort_prefix=cohort_prefix,
            genome_build=genome_build,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_get_call_rate_mt
        }
    }
    call helpers.getHailMTSizes as getCallRateMTSizes {
        input:
        mt_uris=getCallRateMT.call_rate_mt,
        hail_docker=hail_docker
    }

    call inferPlatformPCA {
        input:
            call_rate_mts=getCallRateMT.call_rate_mt,
            hail_docker=hail_docker,
            genome_build=genome_build,
            cohort_set_id=cohort_set_id,
            bucket_id=bucket_id,
            input_size=getCallRateMTSizes.mt_size,
            n_pcs=n_pcs,
            binarization_threshold=binarization_threshold,
            hdbscan_min_cluster_size=hdbscan_min_cluster_size,
            hdbscan_min_samples=hdbscan_min_samples,
            runtime_attr_override=runtime_attr_infer_platform
    }

    output {
        String loadings_ht = inferPlatformPCA.loadings_ht
        String scores_ht =inferPlatformPCA.scores_ht
        File platform_tsv = inferPlatformPCA.platform_tsv
    }
}   

task getCallRateMT {
    input {
        File vcf_file
        File interval_list
        String cohort_prefix
        String genome_build
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    command <<<
    cat <<EOF > get_call_rate_mt.py
    import hail as hl
    import pandas as pd
    import numpy as np
    import os
    import sys
    import datetime
    from gnomad.sample_qc.ancestry import apply_onnx_classification_model, apply_sklearn_classification_model, assign_population_pcs
    from gnomad.utils.filtering import filter_to_adj
    import gnomad.sample_qc.platform 
    import logging
    from typing import List, Optional, Tuple
    from logging import INFO, DEBUG
    logger = logging.getLogger()
    logger.setLevel(INFO)

    vcf_file = sys.argv[1]
    interval_list = sys.argv[2]
    cohort_prefix = sys.argv[3]
    genome_build = sys.argv[4]
    cores = sys.argv[5]
    mem = int(np.floor(float(sys.argv[6])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, reference_genome=genome_build, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'],
                        array_elements_required=False, call_fields=[], find_replace=('nul', '.'))
    intervals = hl.import_locus_intervals(interval_list, reference_genome=genome_build)

    mt = mt.annotate_entries(GT=hl.or_missing(mt.DP>0, mt.GT))
    call_rate_mt = gnomad.sample_qc.platform.compute_callrate_mt(mt, intervals)
    call_rate_mt.write(f"gs://fc-16559487-7c8d-477e-8bce-47462ad86a98/hail/call_rate_mt/{cohort_prefix}_call_rate.mt",
                    overwrite=True)    
    EOF

    python3 get_call_rate_mt.py ~{vcf_file} ~{interval_list} ~{cohort_prefix} ~{genome_build} \
        ~{cpu_cores} ~{memory}
    >>>

    output {
        String call_rate_mt = "gs://fc-16559487-7c8d-477e-8bce-47462ad86a98/hail/call_rate_mt/~{cohort_prefix}_call_rate.mt"
    }
}

task inferPlatformPCA {
    input {
        Array[String] call_rate_mts
        String cohort_set_id
        String hail_docker
        String genome_build
        String bucket_id
        Float input_size
        Float binarization_threshold
        Int n_pcs
        Int hdbscan_min_cluster_size
        Int hdbscan_min_samples
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    cat <<EOF > infer_platform.py
    import hail as hl
    import pandas as pd
    import numpy as np
    import os
    import sys
    import datetime
    from gnomad.sample_qc.ancestry import apply_onnx_classification_model, apply_sklearn_classification_model, assign_population_pcs
    from gnomad.utils.filtering import filter_to_adj
    import gnomad.sample_qc.platform 
    import logging
    from typing import List, Optional, Tuple
    from logging import INFO, DEBUG
    logger = logging.getLogger()
    logger.setLevel(INFO)

    call_rate_mts = sys.argv[1].split(',')
    genome_build = sys.argv[2]
    cohort_set_id = sys.argv[3]
    bucket_id = sys.argv[4]
    n_pcs = int(sys.argv[5])
    hdbscan_min_cluster_size = int(sys.argv[6])
    hdbscan_min_samples = int(sys.argv[7])
    binarization_threshold = float(sys.argv[8])
    cores = sys.argv[9]
    mem = int(np.floor(float(sys.argv[10])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")


    if hdbscan_min_cluster_size==0:
        hdbscan_min_cluster_size = None
    
    def assign_platform_from_pcs(
        platform_pca_scores_ht: hl.Table,
        pc_scores_ann: str = "scores",
        hdbscan_min_cluster_size: Optional[int] = None,
        hdbscan_min_samples: int = 50,  # changed this arbitrarily for now
    ) -> hl.Table:
        """
        Assign platforms using HBDSCAN on the results of call rate PCA.

        :param platform_pca_scores_ht: Input table with the PCA score for each sample
        :param pc_scores_ann: Field containing the scores
        :param hdbscan_min_cluster_size: HDBSCAN `min_cluster_size` parameter. If not specified the smallest of 500 and 0.1*n_samples will be used.
        :param hdbscan_min_samples: HDBSCAN `min_samples` parameter
        :return: A Table with a `qc_platform` annotation containing the platform based on HDBSCAN clustering
        """
        import hdbscan

        logger.info("Assigning platforms based on platform PCA clustering")

        # Read and format data for clustering
        data = platform_pca_scores_ht.to_pandas()
        callrate_data = np.asarray(data[pc_scores_ann].tolist(), dtype=float)
        logger.info("Assigning platforms to %d samples.", len(callrate_data))

        # Cluster data
        if hdbscan_min_cluster_size is None:
            hdbscan_min_cluster_size = min(500, int(0.1 * data.shape[0]))  # edited to int
        
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples
        )
        cluster_labels = clusterer.fit_predict(callrate_data)
        n_clusters = len(set(cluster_labels)) - (
            -1 in cluster_labels
        )  # NOTE: -1 is the label for noisy (un-classifiable) data points
        logger.info("Found %d unique platforms during platform imputation.", n_clusters)

        data["qc_platform"] = cluster_labels
        ht = hl.Table.from_pandas(data, key=[*platform_pca_scores_ht.key])
        ht = ht.annotate(qc_platform="platform_" + hl.str(ht.qc_platform))
        return ht

    for i, uri in enumerate(call_rate_mts):
        if i==0:
            call_rate_mt = hl.read_matrix_table(uri)
        else:
            mt = hl.read_matrix_table(uri)
            call_rate_mt = call_rate_mt.union_cols(mt)
    
    eigenvalues, scores_ht, loadings_ht = gnomad.sample_qc.platform.run_platform_pca(call_rate_mt, n_pcs=n_pcs,
            binarization_threshold=binarization_threshold)
    loadings_ht.write(f"{bucket_id}/hail/infer_platform_pca/{cohort_set_id}_platform_pca_loadings.ht", overwrite=True)
    scores_ht.write(f"{bucket_id}/hail/infer_platform_pca/{cohort_set_id}_platform_pca_scores.ht", overwrite=True)

    platform_ht = assign_platform_from_pcs(scores_ht.annotate(scores=scores_ht.scores[:n_pcs]), 
                                                                    pc_scores_ann='scores', 
                                                    hdbscan_min_cluster_size=None, 
                                                    hdbscan_min_samples=None)
    platform_df = platform_ht.to_pandas()
    platform_df.to_csv(f"{cohort_set_id}_assigned_platforms.tsv", sep='\t', index=False)
    EOF
    python3 infer_platform.py ~{sep=',' call_rate_mts} ~{genome_build} ~{cohort_set_id} ~{bucket_id} \
        ~{n_pcs} ~{hdbscan_min_cluster_size} ~{hdbscan_min_samples} ~{binarization_threshold} ~{cpu_cores} ~{memory}
    >>>

    output {
        String loadings_ht = "~{bucket_id}/hail/infer_platform_pca/~{cohort_set_id}_platform_pca_loadings.ht"
        String scores_ht = "~{bucket_id}/hail/infer_platform_pca/~{cohort_set_id}_platform_pca_scores.ht"
        File platform_tsv = "~{cohort_set_id}_assigned_platforms.tsv"
    }
}