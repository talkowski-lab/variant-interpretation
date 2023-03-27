# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs.wdl"

workflow RdTestVisualization{
    input{
        String prefix
        Array[File] medianfile
        File pedfile
        File sample_batches
        File batch_bincov
        File bed
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_rdtest
    }
        call rdtest{
            input:
                bed=bed,
                medianfile=medianfile,
                pedfile=pedfile,
                sample_batches=sample_batches,
                batch_bincov=batch_bincov,
                prefix=prefix,
                sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                runtime_attr_override = runtime_attr_rdtest
        }
        output{
            File Plots = rdtest.plots
            File allcovfile = rdtest.allcovfile
            File median_file = rdtest.median_file
            File samples_text = rdtest.samples_text
            File test_bed = rdtest.test_bed
        }
}


# Run rdtest
task rdtest {
    input{
        File bed
        File sample_batches # samples, batches
        File batch_bincov # batch, bincov, index
        Array[File] medianfile
        File pedfile
        String prefix
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, sample_batches, batch_bincov, medianfile, pedfile]), "GB")
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
        cat ~{bed} | gunzip | cut -f6 > sample.bed
        cat ~{bed} | gunzip | cut -f1-4 > start.bed
        cat ~{bed} | gunzip | cut -f5 > svtype.bed
        paste start.bed sample.bed svtype.bed > final.bed
        cat final.bed |egrep "DEL|DUP" | sort -k1,1 -k2,2n> test.bed
        cut -f5 test.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        fgrep -wf samples.txt ~{sample_batches} |awk '{print $2}' |sort -u >existing_batches.txt
        fgrep -f existing_batches.txt ~{batch_bincov} > bincovlist.txt
        paste ~{sep=" " medianfile} > medianfile.txt

        i=0
        bedtools merge -i test.bed > test.merged.bed
        while read batch bincov
        do
            let "i=$i+1"
            if [ $i -gt 1 ]
            then
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R test.merged.bed|cut -f4->covfile.$i.bed 
            else
                export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
                tabix -h $bincov -R test.merged.bed>covfile.$i.bed 

        fi
        done<bincovlist.txt

        paste covfile.*.bed |tr ' ' '\t' |bgzip >allcovfile.bed.gz 
        tabix allcovfile.bed.gz
        rm covfile.*.bed
        zcat allcovfile.bed.gz |head -n 1|cut -f 4-|tr '\t' '\n'>samples.txt
        Rscript /opt/RdTest/Rd.R \
            -b test.bed \
            -n ~{prefix} \
            -c allcovfile.bed.gz \
            -m medianfile.txt \
            -f ~{pedfile} \
            -a TRUE \
            -d TRUE \
            -w samples.txt \
            -s 10000000
        mkdir ~{prefix}_rd_plots
        mv *jpg ~{prefix}_rd_plots
        tar -czvf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots/
    >>>
    
    output {
        File plots = "~{prefix}_rd_plots.tar.gz"
        File allcovfile = "allcovfile.bed.gz"
        File median_file = "medianfile.txt"
        File test_bed = "test.bed"
        File samples_text = "samples.txt"
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