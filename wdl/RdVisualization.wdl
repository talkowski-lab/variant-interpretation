# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs.wdl"

workflow RdTestVisualization{
    input{
        String prefix
        File? fam_ids
        Array[File] medianfile
        File pedfile
        File sample_batches
        File batch_bincov
        File bed
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_rdtest
    }

    if (defined(fam_ids)) {
        File fam_ids_ = select_first([fam_ids])
        Array[String] family_ids = transpose(read_tsv(fam_ids_))[0]
    }

    if (!(defined(fam_ids))) {
        call generate_families{
            input:
                bed = bed,
                ped_file = pedfile,
                sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
                runtime_attr_override = runtime_attr_rdtest
        }
    }

    scatter (family in select_first([family_ids, generate_families.families])){
        call rdtest{
            input:
                bed=bed,
                family = family,
                ped_file = pedfile,
                medianfile=medianfile,
                pedfile=pedfile,
                sample_batches=sample_batches,
                batch_bincov=batch_bincov,
                prefix=prefix,
                sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                runtime_attr_override = runtime_attr_rdtest
        }
    }

    call integrate_rd_plots{
        input:
            rd_tar = rdtest.plots,
            prefix = prefix, 
            sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
            runtime_attr_override = runtime_attr_rdtest
    }
        output{
            File Plots = integrate_rd_plots.plot_tar
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
        String family
        File ped_file
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
        cat ~{ped_file} | grep -w ~{family} | cut -f2 | sort -u > samples_in_family.txt
        cat ~{bed} | gunzip | cut -f1-6 | grep -w -f samples_in_family.txt > per_family_bed.bed
        cat per_family_bed.bed | gunzip | cut -f6 > sample.bed
        cat per_family_bed.bed | gunzip | cut -f1-4 > start.bed
        cat per_family_bed.bed | gunzip | cut -f5 > svtype.bed
        paste start.bed sample.bed svtype.bed > final.bed
        cat final.bed |egrep "DEL|DUP" | sort -k1,1 -k2,2n> test.bed
        cut -f5 test.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        cat ~{ped_file} | grep -w -f samples.txt | cut -f1 | sort -u > families.txt
        cat ~{ped_file} | grep -w -f families.txt | cut -f2 | sort -u > all_samples.txt
        fgrep -wf all_samples.txt ~{sample_batches} |awk '{print $2}' |sort -u >existing_batches.txt
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
        mkdir rd_plots
        zcat allcovfile.bed.gz |head -n 1|cut -f 4-|tr '\t' '\n'>samples.txt
        Rscript /opt/RdTest/Rd.R \
            -b test.bed \
            -n ~{family} \
            -c allcovfile.bed.gz \
            -m medianfile.txt \
            -f ~{pedfile} \
            -a TRUE \
            -d TRUE \
            -w samples.txt \
            -s 10000000
        tar -czvf ~{family}_rd_plots.tar.gz rd_plots
    >>>
    
    output {
        File plots = "~{family}_rd_plots.tar.gz"
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

task generate_families{
    input {
        File bed
        File ped_file
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, ped_file]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    command <<<
        set -euo pipefail
        cat ~{bed} | gunzip | tail -n+2 | cut -f6 | tr ',' '\n' | sort -u > samples.txt #must have header line
        grep -w -f samples.txt ~{ped_file} | cut -f1 | sort -u  > families.txt
        >>>

    output{
        Array[String] families = read_lines("families.txt")
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: sv_pipeline_rdtest_docker
    }

}

task integrate_rd_plots{
    input {
        Array[File] rd_tar
        String prefix
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(rd_tar, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    command <<<
        mkdir ~{prefix}_rd_plots
        while read file; do
            tar -zxf ${file}
            mv pe_rd_plots/*  ~{prefix}_rd_plots/
        done < ~{write_lines(rd_tar)};
        tar -czf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots
    >>>

    output{
        File plot_tar = "~{prefix}_igv_plots.tar.gz"
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