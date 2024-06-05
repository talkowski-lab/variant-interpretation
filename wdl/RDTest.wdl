# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs2.wdl"

workflow RdTest{
    input{
        String prefix
        Array[File] medianfile
        File sample_batches
        File outlier_samples
        File batch_bincov
        File bed
        Int rd_window
        String sv_pipeline_rdtest_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_rdtest
        RuntimeAttr? runtime_attr_create_bed
    }

    call rdtest{
        input:
            bed = bed,
            medianfile = medianfile,
            sample_batches = sample_batches,
            outlier_samples = outlier_samples,
            batch_bincov = batch_bincov,
            prefix = prefix,
            rd_window = rd_window,
            sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
            runtime_attr_override = runtime_attr_rdtest
    }

    output{
        File Plots = rdtest.plots
    }
}

# Run rdtest
task rdtest {
    input{
        File bed
        File sample_batches # samples, batches
        File batch_bincov # batch, bincov, index
        Array[File] medianfile
        File outlier_samples
        Int rd_window
        String prefix
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, sample_batches, batch_bincov, ped_file]), "GB")
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
        cat ~{bed} |egrep "DEL|DUP" | sort -k1,1 -k2,2n> input.bed
        cut -f5 input.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        fgrep -wf samples.txt ~{sample_batches} |awk '{print $2}' |sort -u > existing_batches.txt
        grep -w -f existing_batches.txt ~{batch_bincov} > bincovlist.txt
        paste ~{sep=" " medianfile} > medianfile.txt

        i=0
        bedtools merge -i input.bed > input.merged.bed
        while read batch bincov index
        do
            let "i=$i+1"
            if [ $i -gt 1 ]
            then
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R input.merged.bed | cut -f4- > covfile.$i.bed
            else
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R input.merged.bed > covfile.$i.bed

        fi
        done<bincovlist.txt

        paste covfile.*.bed | tr ' ' '\t' | bgzip > allcovfile.bed.gz
        tabix allcovfile.bed.gz
        rm covfile.*.bed
        zcat allcovfile.bed.gz | head -n 1 | cut -f 4- | tr '\t' '\n' > all_samples.txt

        ##Remove outliers from RD plot except samples with the variant
        grep -vf samples.txt ~{outlier_samples} > outliers_keep.txt
        grep -vf outliers_keep.txt allcovfile.bed.gz > samples_noOutliers.txt

        ##Run RD test script
        Rscript /opt/RdTest/Rd.R \
            -b input.bed \
            -n ~{prefix} \
            -c allcovfile.bed.gz \
            -m medianfile.txt \
            -a TRUE \
            -d TRUE \
            -w samples_noOutliers.txt \
            -s ~{rd_window}

        mkdir rd_plots
        mv *jpg rd_plots
        tar -czvf rd_plots.tar.gz rd_plots/
    >>>

    output {
        File plots = "rd_plots.tar.gz"
        File allcovfile = "allcovfile.bed.gz"
        File median_file = "medianfile.txt"
        File test_bed = "input.bed"
        File samples_text = "samples_noOutliers.txt"
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