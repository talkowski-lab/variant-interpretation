# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs2.wdl"

workflow RdTestVisualization{
    input{
        String prefix
        File? fam_ids
        File batch_medianfile
        File ped_file
        File sample_batches
        File outlier_samples
        File batch_bincov
        File bed
        File? regeno
        String sv_pipeline_rdtest_docker
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_rdtest
        RuntimeAttr? runtime_attr_create_bed
    }

    if (defined(fam_ids)) {
        File fam_ids_ = select_first([fam_ids])
        Array[String] family_ids = transpose(read_tsv(fam_ids_))[0]
    }
    #TO FIX: Currently, if defined families, rd fails is a family does not have DEL/DUP
    if (!(defined(fam_ids))) {
        call generate_families{
            input:
                bed = bed,
                ped_file = ped_file,
                sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
                runtime_attr_override = runtime_attr_rdtest
        }
    }

    scatter (family in select_first([family_ids, generate_families.families])){
        call generatePerFamilyBed{
            input:
                bed=bed,
                family = family,
                batch_medianfile = batch_medianfile,
                sample_batches = sample_batches,
                ped_file = ped_file,
                variant_interpretation_docker = variant_interpretation_docker,
                runtime_attr_override = runtime_attr_create_bed
        }

        if (defined(regeno)){
            File regeno_ = select_first([regeno])
            call rdtest_regeno{
                input:
                    bed=generatePerFamilyBed.bed_file,
                    family = family,
                    ped_file = ped_file,
                    medianfile = generatePerFamilyBed.medianfile,
                    sample_batches=sample_batches,
                    outlier_samples=outlier_samples,
                    regeno = regeno_,
                    batch_bincov=batch_bincov,
                    prefix=prefix,
                    sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                    runtime_attr_override = runtime_attr_rdtest
                }
        }

        if (!defined(regeno)){
            call rdtest{
                input:
                    bed=generatePerFamilyBed.bed_file,
                    family = family,
                    ped_file = ped_file,
                    medianfile = generatePerFamilyBed.medianfile,
                    sample_batches=sample_batches,
                    outlier_samples=outlier_samples,
                    batch_bincov=batch_bincov,
                    prefix=prefix,
                    sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                    runtime_attr_override = runtime_attr_rdtest
                }
        }
    }

    call integrate_rd_plots{
        input:
#            rd_tar = select_first([rdtest_regeno.plots, rdtest.plots]),
            rd_tar = if (defined (regeno)) then select_all(rdtest_regeno.plots) else select_all(rdtest.plots),
            prefix = prefix, 
            sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
            runtime_attr_override = runtime_attr_rdtest
    }

    output{
        File Plots = integrate_rd_plots.plot_tar
        File median = rdtest.median_file
        File? median_geno = rdtest_regeno.rd_median_geno
    }
}

task generatePerFamilyBed {
    input{
        File bed
        String family
        File ped_file
        File batch_medianfile
        File sample_batches
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([bed, ped_file, batch_medianfile, sample_batches]), "GB")
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
        cat ~{bed} | gunzip | cut -f1-6 | bgzip -c > first_six.bed.gz
        Rscript /src/variant-interpretation/scripts/reformatSingleSampleBed.R first_six.bed.gz
        cat single_sample.bed | grep -w -f samples_in_family.txt > per_family_bed.bed
        cat per_family_bed.bed | cut -f6 > sample.bed
        cat per_family_bed.bed | cut -f1-4 > start.bed
        cat per_family_bed.bed | cut -f5 > svtype.bed
        paste start.bed sample.bed svtype.bed > final.bed
        grep -w -f samples_in_family.txt ~{sample_batches} |awk '{print $2}' |sort -u >existing_batches.txt
        grep -w -f existing_batches.txt ~{batch_medianfile} | cut -f2 > medianfile.txt
    >>>
    
    output {
        File bed_file = "final.bed"
        Array[File] medianfile = read_lines("medianfile.txt")
    }
    
    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: variant_interpretation_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        File outlier_samples
        Array[File] medianfile
        String prefix
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, sample_batches, batch_bincov, medianfile, ped_file]), "GB")
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
        cat ~{bed} |egrep "DEL|DUP" | sort -k1,1 -k2,2n> test.bed
        cut -f5 test.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        cat ~{ped_file} | grep -w -f samples.txt | cut -f1 | sort -u > families.txt
        cat ~{ped_file} | grep -w -f families.txt | cut -f2 | sort -u > all_samples.txt
        fgrep -wf all_samples.txt ~{sample_batches} |awk '{print $2}' |sort -u >existing_batches.txt
        grep -w -f existing_batches.txt ~{batch_bincov} > bincovlist.txt
        paste ~{sep=" " medianfile} > medianfile.txt

        i=0
        bedtools merge -i test.bed > test.merged.bed
        while read batch bincov index
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

        ##Pass only subset ped file
        head -n+1 ~{ped_file} > subset_families.ped
        grep -wf families.txt ~{ped_file} >> subset_families.ped
        cut -f2 subset_families.ped | sort -u > subset_samples.txt

        ##Remove outliers from RD plot except samples from this family
        grep -vf subset_samples.txt ~{outlier_samples} > outliers_keepFams.txt
        grep -vf outliers_keepFams.txt samples.txt > samples_noOutliers.txt

        ##Run RD test script
        Rscript /opt/RdTest/Rd.R \
            -b test.bed \
            -n ~{prefix} \
            -c allcovfile.bed.gz \
            -m medianfile.txt \
            -f subset_families.ped \
            -a TRUE \
            -d TRUE \
            -w samples_noOutliers.txt \
            -s 10000000

        mkdir rd_plots
        mv *jpg rd_plots
        tar -czvf rd_plots.tar.gz rd_plots/
    >>>

    output {
        File plots = "rd_plots.tar.gz"
        File allcovfile = "allcovfile.bed.gz"
        File median_file = "medianfile.txt"
        File test_bed = "test.bed"
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

# Run rdtest
task rdtest_regeno {
    input{
        File bed
        String family
        File ped_file
        File sample_batches # samples, batches
        File batch_bincov # batch, bincov, index
        File outlier_samples
        File regeno
        Array[File] medianfile
        String prefix
        String sv_pipeline_rdtest_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(select_all([bed, sample_batches, batch_bincov, medianfile, ped_file, regeno, outlier_samples]), "GB")
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
        cat ~{bed} |egrep "DEL|DUP" | sort -k1,1 -k2,2n> test.bed
        cut -f5 test.bed |sed 's/\,/\n/g'|sort -u > samples.txt
        cat ~{ped_file} | grep -w -f samples.txt | cut -f1 | sort -u > families.txt
        cat ~{ped_file} | grep -w -f families.txt | cut -f2 | sort -u > all_samples.txt
        fgrep -wf all_samples.txt ~{sample_batches} |awk '{print $2}' |sort -u >existing_batches.txt
        grep -w -f existing_batches.txt ~{batch_bincov} > bincovlist.txt
        paste ~{sep=" " medianfile} > medianfile.txt

        i=0
        bedtools merge -i test.bed > test.merged.bed
        while read batch bincov index
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

        ##Pass only subset ped file
        head -n+1 ~{ped_file} > subset_families.ped
        grep -wf families.txt ~{ped_file} >> subset_families.ped
        cut -f2 subset_families.ped | sort -u > subset_samples.txt

        ##Remove outliers from RD plot except samples from this family
        grep -vf subset_samples.txt ~{outlier_samples} > outliers_keepFams.txt
        grep -vf outliers_keepFams.txt samples.txt > samples_noOutliers.txt

        ##Run RD test script
        Rscript /opt/RdTest/Rd.R \
            -b test.bed \
            -n ~{prefix} \
            -c allcovfile.bed.gz \
            -m medianfile.txt \
            -f subset_families.ped \
            -a TRUE \
            -d TRUE \
            -w samples_noOutliers.txt \
            -s 10000000 \
            -g TRUE
            -r ~{regeno}

        mkdir rd_plots
        mv *jpg rd_plots
        tar -czvf rd_plots.tar.gz rd_plots/
    >>>
    
    output {
        File plots = "rd_plots.tar.gz"
        File allcovfile = "allcovfile.bed.gz"
        File median_file = "medianfile.txt"
        File test_bed = "test.bed"
        File samples_text = "samples_noOutliers.txt"
        File rd_median_geno = "~{prefix}.median_geno"
        File rd_denovo = "~{prefix}.denovo"
        File rd_geno = "~{prefix}.geno"
        File rd_gq = "~{prefix}.gq"
        File rd_vargq = "~{prefix}.vargq"
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
        cat ~{bed} | gunzip | tail -n+2 | cut -f 1-6 | grep 'DEL\|DUP' | cut -f6 | tr ',' '\n' | sort -u > samples.txt #must have header line
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
            mv rd_plots/*  ~{prefix}_rd_plots/
        done < ~{write_lines(rd_tar)};
        
        tar -czf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots
    >>>

    output{
        File plot_tar = "~{prefix}_rd_plots.tar.gz"
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