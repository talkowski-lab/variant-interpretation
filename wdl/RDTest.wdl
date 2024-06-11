# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs2.wdl"

workflow RdTest{
    input{
        String prefix
        File medianfile
        File sample_batches
        File outlier_samples
        File batch_bincov
        File bed
        File ped_file
        Int rd_window
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_rdtest
    }

    Array[String] variants = transpose(read_tsv(bed))[3]

    scatter(var in variants) {
        call rdtest{
            input:
                bed = bed,
                variant = var,
                medianfile = medianfile,
                sample_batches = sample_batches,
                outlier_samples = outlier_samples,
                batch_bincov = batch_bincov,
                prefix = prefix,
                ped_file = ped_file,
                rd_window = rd_window,
                sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
                runtime_attr_override = runtime_attr_rdtest
        }
    }

    output{
        Array[File] Plots = rdtest.plots
    }
}

# Run rdtest
task rdtest {
    input{
        File bed
        File sample_batches # samples, batches
        File batch_bincov # batch, bincov, index
        File medianfile
        File ped_file
        File outlier_samples
        String variant
        Int rd_window
        String prefix
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, sample_batches, batch_bincov]), "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(base_disk_gb + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -ex
        cat ~{bed} | grep -w ~{variant} | sed -e 's/CPX/DEL/g' | sort -k1,1 -k2,2n> input.bed
        cut -f5 input.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        fgrep -wf samples.txt ~{sample_batches} |awk '{print $2}' |sort -u > existing_batches.txt
        grep -w -f existing_batches.txt ~{batch_bincov} > bincovlist.txt

        i=0

        while read batch bincov index
        do
            let "i=$i+1"
            if [ $i -gt 1 ]
            then
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R input.bed | cut -f4- > covfile.$i.bed
            else
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R input.bed > covfile.$i.bed

        fi
        done<bincovlist.txt

        paste covfile.*.bed | tr ' ' '\t' | sort -k1,1 -k2,2n | bgzip > allcovfile.bed.gz
        tabix allcovfile.bed.gz
        rm covfile.*.bed
        zcat allcovfile.bed.gz | head -n 1 | cut -f 4- | tr '\t' '\n' > all_samples.txt

        ##Remove outliers from RD plot except samples with the variant
        grep -vf samples.txt ~{outlier_samples} > outliers_keep.txt
        grep -vf outliers_keep.txt all_samples.txt > samples_noOutliers.txt

        ##Make output directory
        mkdir ~{variant}

        ##Run RD test script
        Rscript /opt/RdTest/RdTest.R \
            -b input.bed \
            -n ~{prefix} \
            -c allcovfile.bed.gz \
            -m ~{medianfile} \
            -w samples_noOutliers.txt \
            -s ~{rd_window} \
            -f ~{ped_file} \
            -p TRUE \
            -o ~{variant}

        tar -czvf ~{variant}.tar.gz ~{variant}
    >>>

    output {
        File plots = "~{variant}.tar.gz"
#        File allcovfile = "allcovfile.bed.gz"
#        File test_bed = "input.bed"
#        File samples_text = "samples_noOutliers.txt"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_rdtest_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}