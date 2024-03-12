version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

task getHailMTSize {
    input {
        String mt_uri
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
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
        gsutil du -sh ~{mt_uri} | cut -f1 -d ' ' > mt_size.txt
    >>>

    output {
        Float mt_size = read_lines('mt_size.txt')[0]
    }
}

task getHailMTSizes {
    input {
        Array[String] mt_uris
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
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
        touch mt_size.txt
        for mt_uri in $(cat ~{write_lines(mt_uris)});
        do
            gsutil -m du -sh $mt_uri | cut -f1 -d ' ' >> mt_size.txt;
        done

        cat <<EOF > get_sum.py
        import pandas as pd
        sizes = pd.read_csv('mt_size.txt', header=None)[0]
        pd.Series([sizes.sum()]).to_csv('tot_size.txt', index=False, header=None)
        EOF

        python3 get_sum.py
    >>>

    output {
        Float mt_size = read_lines('tot_size.txt')[0]
    }
}

task filterIntervalsToMT {
    input {
        File bed_file
        Float input_size
        String mt_uri
        String bucket_id        
        String hail_docker
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
    cat <<EOF > filter_intervals.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uri = sys.argv[1]
    bed_file = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    bucket_id = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.read_matrix_table(mt_uri)
    intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
    mt_filt = mt.filter_rows(hl.is_defined(intervals[mt.locus]))

    filename = f"{bucket_id}/hail/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{os.path.basename(mt_uri).split('.mt')[0]}_{os.path.basename(bed_file).split('.bed')[0]}.mt"
    mt_filt.write(filename)
    pd.Series([filename]).to_csv('mt_uri.txt', index=False, header=None)
    EOF
    
    python3 filter_intervals.py ~{mt_uri} ~{bed_file} ~{cpu_cores} ~{memory} ~{bucket_id}
    >>>

    output {
        String mt_filt = read_lines('mt_uri.txt')[0]
    }
}

task filterIntervalsToVCF {
    input {
        File bed_file
        Float input_size
        String mt_uri
        String hail_docker
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
    cat <<EOF > filter_intervals.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uri = sys.argv[1]
    bed_file = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.read_matrix_table(mt_uri)
    intervals = hl.import_bed(bed_file, reference_genome='GRCh38')
    mt_filt = mt.filter_rows(hl.is_defined(intervals[mt.locus]))

    filename = f"{os.path.basename(mt_uri).split('.mt')[0]}_{os.path.basename(bed_file).split('.bed')[0]}.vcf.bgz"
    hl.export_vcf(mt_filt, filename)
    EOF
    
    python3 filter_intervals.py ~{mt_uri} ~{bed_file} ~{cpu_cores} ~{memory}
    >>>

    output {
        File vcf_filt = "~{basename(mt_uri, '.mt')}_~{basename(bed_file, '.bed')}.vcf.bgz"
    }
}

task subsetVCFSamplesHail {
    input {
        File vcf_file
        File samples_file  # .txt extension  
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
    cat <<EOF > subset_samples.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    samples_file = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_file, reference_genome = 'GRCh38', array_elements_required=False, force_bgz=True)
    samples = pd.read_csv(samples_file, header=None)[0].tolist()

    mt_filt = mt.filter_cols(hl.array(samples).contains(mt.s))
    hl.export_vcf(mt_filt, os.path.basename(samples_file).split('.txt')[0]+'.vcf.bgz')
    EOF

    python3 subset_samples.py ~{vcf_file} ~{samples_file} ~{cpu_cores} ~{memory}
    >>>

    output {
        File vcf_subset = basename(samples_file, '.txt') + '.vcf.bgz'
    }
}

task mergeMTs {
    input {
        Array[String] mt_uris
        String cohort_prefix
        String bucket_id
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
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
    cat <<EOF > merge_mts.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    mt_uris = sys.argv[1].split(',')
    merged_filename = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))
    bucket_id = sys.argv[5]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    tot_mt = len(mt_uris)
    for i, mt_uri in enumerate(mt_uris):
        if (((i+1)%10)==0):
            print(f"Merging MT {i+1}/{tot_mt}...")
        if i==0:
            mt = hl.read_matrix_table(mt_uri)
        else:
            mt2 = hl.read_matrix_table(mt_uri)
            mt = mt.union_rows(mt2)
    filename = f"{bucket_id}/hail/merged_mt/{str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M'))}/{merged_filename}.mt"
    mt.write(filename, overwrite=True)
    pd.Series([filename]).to_csv('mt_uri.txt', index=False, header=None)
    EOF

    python3 merge_mts.py ~{sep=',' mt_uris} ~{cohort_prefix}_merged ~{cpu_cores} ~{memory} ~{bucket_id}
    >>>

    output {
        String merged_mt = read_lines('mt_uri.txt')[0]
    }
}

task mergeHTs {
    input {
        Array[String] ht_uris
        String merged_filename
        String hail_docker
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
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
    cat <<EOF > merge_hts.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    ht_uris = sys.argv[1].split(',')
    merged_filename = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{mem}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{mem}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    tot_ht = len(ht_uris)
    merged_df = pd.DataFrame()

    for i, uri in enumerate(ht_uris):
        print(f"{i+1}/{tot_ht}")
        ht = hl.read_table(uri)
        merged_df = pd.concat([merged_df, ht.to_pandas()])

    merged_df.to_csv(f"{merged_filename}.tsv", sep='\t', index=False)
    EOF

    python3 merge_hts.py ~{sep=',' ht_uris} ~{merged_filename} ~{cpu_cores} ~{memory}
    >>>

    output {
        File merged_tsv = merged_filename + '.tsv'
    }
}