version 1.0
    
import "Structs.wdl"
import "SplitVcf.wdl" as splitVCF
import "ReformatRawFiles.wdl" as reformatRaw
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "DeNovoSVsScatter.wdl" as runDeNovo

workflow GenomicDisorders {

    input {
        File ped_input
        File vcf_file
        Array[String] contigs
        File genomic_disorder_input
        File genomic_disorder_input_ref
        File batch_raw_file
        File batch_depth_raw_file
        File sample_batches
        File batch_bincov_index
        String prefix
        File? fam_ids
        String variant_interpretation_docker
        String genomic_disorders_docker
        RuntimeAttr? runtime_attr_get_batched_files
        RuntimeAttr? runtime_attr_clean_ped
        RuntimeAttr? runtime_attr_raw_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_merge_bed
        RuntimeAttr? runtime_attr_raw_divide_by_chrom
        RuntimeAttr? runtime_attr_raw_reformat_bed
        RuntimeAttr? runtime_attr_gd
        RuntimeAttr? runtime_attr_merge_gd
        RuntimeAttr? runtime_attr_reformat_vcf
        RuntimeAttr? runtime_attr_vcf_overlap
    }

    #if the fam_ids input is given, subset all other input files to only include the necessary batches
    if (defined(fam_ids)){
        File fam_ids_ = select_first([fam_ids])
        call getBatchedFiles{
            input:
                batch_raw_file = batch_raw_file,
                batch_depth_raw_file = batch_depth_raw_file,
                ped_input = ped_input,
                fam_ids = fam_ids_,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_get_batched_files
        }
    }
    #Makes a ped file of singletons, duos, and trios for input into the de novo GD filtering
    call cleanPed{
        input:
            ped_input = ped_input,
            vcf_input = vcf_file,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_clean_ped
    }

    #Splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call reformatRaw.reformatRawFiles as reformatRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([getBatchedFiles.batch_raw_files_list, batch_raw_file]),
            ped_input = cleanPed.cleaned_ped,
            depth = false,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }

    #Splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call reformatRaw.reformatRawFiles as reformatDepthRawFiles {
        input:
            contigs = contigs,
            raw_files_list = select_first([getBatchedFiles.batch_depth_raw_files_list, batch_depth_raw_file]),
            ped_input = cleanPed.cleaned_ped,
            depth = true,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_vcf_to_bed = runtime_attr_raw_vcf_to_bed,
            runtime_attr_merge_bed = runtime_attr_raw_merge_bed,
            runtime_attr_divide_by_chrom = runtime_attr_raw_divide_by_chrom,
            runtime_attr_reformat_bed = runtime_attr_raw_reformat_bed
    }
    
    scatter (i in range(length(contigs))){
        #Generates a list of genomic disorder regions in the vcf input as well as in the depth raw files
        call getGDraw{
            input:
                genomic_disorder_input=genomic_disorder_input,
                ped = ped_input,
#                vcf_file = select_first([splitVCF.split_vcf, vcf_file]),
                depth_raw_file_proband = reformatDepthRawFiles.reformatted_proband_raw_files[i],
                depth_raw_file_parents = reformatDepthRawFiles.reformatted_parents_raw_files[i],
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_gd
        }

        call getDenovoGDraw{
            input:
                ped = ped_input,
                gd_proband_calls = getGDraw.gd_proband_calls[i],
                gd_parent_calls = getGDraw.gd_parent_calls[i],
                gd_coverage_proband = getGDraw.gd_coverage_proband[i],
                gd_coverage_parents = getGDraw.gd_coverage_parents[i],
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_gd
        }
    }
    #Merges the genomic disorder region output from each chromosome to compile a list of genomic disorder regions
    call mergeGenomicDisorders{
        input:
            gd_bed_to_merge=getGDraw.gd_output_from_depth_raw_files,
            gd_denovo_bed_to_merge=getGDraw.gd_output_from_depth_raw_files_denovo,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_gd
    }

    call reformatVCF{
        input:
            vcf = vcf_file,
            prefix = prefix,
            docker_path = genomic_disorders_docker,
            runtime_attr_override = runtime_attr_reformat_vcf
    }

    call getGDvcf{
        input:
          bed = reformatVCF.out_ref_bed,
          genomic_disorders = genomic_disorder_input_ref,
          prefix = prefix,
          docker_path = genomic_disorders_docker,
          runtime_attr_override = runtime_attr_vcf_overlap
    }

    output {
        File cleaned_ped = cleanPed.cleaned_ped
        File gd_depth = mergeGenomicDisorders.gd_output_from_depth
        File gd_depth_denovo = mergeGenomicDisorders.gd_denovo_output_from_depth
        File vcf_in_gds = getGDvcf.out_bed
    }
}

task getBatchedFiles{
    input{
        File batch_raw_file
        File batch_depth_raw_file
        File fam_ids
        File ped_input
        File sample_batches
        File batch_bincov_index
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([batch_raw_file, batch_depth_raw_file, ped_input, sample_batches, batch_bincov_index, fam_ids]), "GB")
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

    output{
        File batch_raw_files_list = "batch_raw_files_list.txt"
        File batch_depth_raw_files_list = "batch_depth_raw_files_list.txt"
        File batch_bincov_index_subset = "batch_bincov_index.txt"
        File samples = "samples.txt"
    }

    command <<<
        set -euo pipefail
        grep -w -f ${fam_ids} ${ped_input} | cut -f2 | sort -u > samples.txt
        grep -w -f samples.txt ${sample_batches} | cut -f2 | sort -u > batches.txt
        grep -w -f batches.txt ${batch_bincov_index} > batch_bincov_index.txt
        grep -w -f batches.txt ${batch_raw_file} > batch_raw_files_list.txt
        grep -w -f batches.txt ${batch_depth_raw_file} > batch_depth_raw_files_list.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task cleanPed{
    input{
        File ped_input
        File vcf_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float ped_size = size(ped_input, "GB")
    Float vcf_size = size(vcf_input, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + vcf_size + ped_size * 1.5),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File cleaned_ped = "subset_ped.txt"
    }

    command <<<
        set -euo pipefail

        ##Get samples in VCF file
        bcftools query -l ~{vcf_input} > samples_in_vcf.txt

        ##Keep only samples in the VCF file
        grep ^FamID ~{ped_input} > subset_ped.txt
        grep -wf samples_in_vcf.txt ~{ped_input} >> subset_ped.txt

        Rscript /src/variant-interpretation/scripts/cleanPed.R subset_ped.txt

#        cut -f2 cleaned_ped.txt | tail -n+2 > all_samples.txt
#        grep -w -v -f samples_to_include_in_ped.txt all_samples.txt > excluded_samples.txt
#        grep -w -f excluded_samples.txt cleaned_ped.txt | cut -f1 | sort -u > excluded_families.txt
#        grep -w -v -f excluded_families.txt cleaned_ped.txt > subset_cleaned_ped.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task getGDraw{
    input{
        File ped
        File depth_raw_file_proband
        File depth_raw_file_parents
        String chromosome
        File genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([ped, genomic_disorder_input, depth_raw_file_parents, depth_raw_file_proband]), "GB")
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

    output{
        File gd_output_from_depth_raw_files = "~{chromosome}.gd.variants.in.depth.raw.files.txt.gz"
        File gd_proband_calls = "~{chromosome}.proband.GD.calls.txt"
        File gd_parent_calls = "~{chromosome}.parent.GD.calls.txt"
        File gd_coverage_proband = "~{chromosome}.kept.coverage.proband.txt"
        File gd_coverage_parents = "~{chromosome}.kept.coverage.parents.txt"
    }

    command <<<
        set -euxo pipefail

        Rscript /src/variant-interpretation/scripts/create_per_sample_bed.R ~{genomic_disorder_input} unsorted.gd.per.sample.txt unsorted.gd.per.family.txt ~{ped} ~{chromosome}
        sort -k1,1 -k2,2n unsorted.gd.per.sample.txt > gd.per.sample.txt
        sort -k1,1 -k2,2n unsorted.gd.per.family.txt > gd.per.family.txt
        cat ~{depth_raw_file_parents} | gunzip | sort -k1,1 -k2,2n | bgzip -c > sorted.depth.parents.bed.gz
        cat ~{depth_raw_file_proband} | gunzip | sort -k1,1 -k2,2n | bgzip -c > sorted.depth.proband.bed.gz

        echo "Done with R script"

        bedtools intersect -wa -wb -f 0.3 -r -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.proband.txt
        bedtools intersect -wa -wb -f 0.3 -r -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.parents.txt

        echo "done with intersect in depth variants"

        # We should have the same format of bedtools intersect in order to concatanate the files in downstream.
        # The "bedtools coverage" below collect the list of disorders that has hits in depth file limited to coverage fraction.
        # The output file doesn't include the list of sample or family calls which cover 30% of disorder region.
        # Thus, "bedtools intersect" is required to recover the calls.
        bedtools coverage -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz |awk '{if ($NF>=0.30) print }' > ~{chromosome}.coverage.parents.list.txt
        bedtools intersect -wa -wb -a ~{chromosome}.coverage.parents.list.txt -b sorted.depth.parents.bed.gz | sort | uniq | awk -v OFS="\t" '{print $1,$2,$3,$4,$9,$10,$11,$12,$13}' > ~{chromosome}.coverage.parents.txt

        bedtools coverage -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz |awk '{if ($NF>=0.30) print }' > ~{chromosome}.coverage.proband.list.txt
        bedtools intersect -wa -wb -a ~{chromosome}.coverage.proband.list.txt -b sorted.depth.proband.bed.gz | sort | uniq | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$10,$11,$12,$13,$14}' > ~{chromosome}.coverage.proband.txt

        # IMPORTANT "bedtools coverage" throws an error when sorted.depth.proband.bed.gz contains a sample, but the query file "gd.per.sample.txt" (so the ped file) doesn't contain the corresponding sample.
        # This error fix should be done in prior steps.

        bedtools intersect -wa -wb -f 0.3 -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.proband.no.r.txt
        bedtools intersect -wa -wb -f 0.3 -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.parents.no.r.txt

        # Remove the calls that do not reciprocally overlap with the 0.3 fraction of GD site.
        bedtools intersect -v -wb -b ~{chromosome}.gd.variants.in.depth.raw.file.proband.no.r.txt -a ~{chromosome}.coverage.proband.txt > ~{chromosome}.kept.coverage.proband.txt
        bedtools intersect -v -wb -b ~{chromosome}.gd.variants.in.depth.raw.file.parents.no.r.txt -a ~{chromosome}.coverage.parents.txt > ~{chromosome}.kept.coverage.parents.txt
        
        #concatanate proband calls on GD site with desired overlap
        cat ~{chromosome}.gd.variants.in.depth.raw.file.proband.txt ~{chromosome}.kept.coverage.proband.txt | awk -v OFS="\t" '{print $5,$7,$8,$9,$10,$1,$2,$3,$4}' > ~{chromosome}.proband.GD.calls.txt

        #concatanate parent calls on GD site with desired overlap
        cat ~{chromosome}.gd.variants.in.depth.raw.file.parents.txt ~{chromosome}.kept.coverage.parents.txt | awk -v OFS="\t" '{print $5,$6,$7,$8,$9,$1,$2,$3,$4}' > ~{chromosome}.parent.GD.calls.txt

        #to get all GD calls of probands+parents
        cat ~{chromosome}.proband.GD.calls.txt ~{chromosome}.parent.GD.calls.txt > ~{chromosome}.parent_and_proband.GD.calls.txt

        #format the output files: remove sample/family tag from chr in final bed file
        cat ~{chromosome}.parent_and_proband.GD.calls.txt |\
            awk -v OFS="\t" '{sub(/_.*/, "", $1); print}' |\
            awk -v OFS="\t" '{sub(/_.*/, "", $6); print}' | sort | uniq > ~{chromosome}.gd.variants.in.depth.raw.files.txt

        bgzip ~{chromosome}.gd.variants.in.depth.raw.files.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task getDenovoGDraw{
    input{
        File ped
        File gd_proband_calls
        File gd_parent_calls
        File gd_coverage_proband
        File gd_coverage_parents
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([ped, gd_proband_calls, gd_parent_calls]), "GB")
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

    output{
        File gd_output_from_depth_raw_files_denovo ="~{chromosome}.gd.variants.in.depth.raw.files.de.novo.txt.gz"
    }

    command <<<
        set -euxo pipefail

        #concatanate proband calls on GD site with desired overlap
        cat ~{gd_proband_calls} ~{gd_coverage_proband} | awk -v OFS="\t" '{print $5,$7,$8,$9,$10,$1,$2,$3,$4}' > ~{chromosome}.proband.GD.calls.txt

        #concatanate parent calls on GD site with desired overlap
        cat ~{gd_parent_calls} ~{gd_coverage_parents} | awk -v OFS="\t" '{print $5,$6,$7,$8,$9,$1,$2,$3,$4}' > ~{chromosome}.parent.GD.calls.txt

        #get de novo only calls
        bedtools intersect -v -wa  -f 0.3 -a ~{chromosome}.proband.GD.calls.txt -b ~{chromosome}.parent.GD.calls.txt > ~{chromosome}.proband.GD.calls.de_novo.txt

        #to get all GD calls of probands+parents
        cat ~{chromosome}.proband.GD.calls.txt ~{chromosome}.parent.GD.calls.txt > ~{chromosome}.parent_and_proband.GD.calls.txt

        cat ~{chromosome}.proband.GD.calls.de_novo.txt |\
            awk -v OFS="\t" '{sub(/_.*/, "", $1); print}' |\
            awk -v OFS="\t" '{sub(/_.*/, "", $6); print}' | sort | uniq > ~{chromosome}.gd.variants.in.depth.raw.files.de.novo.txt

        bgzip ~{chromosome}.gd.variants.in.depth.raw.files.de.novo.txt
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}


#task subsetVcf{
#    input{
#        File vcf_file
#        String chromosome
#        String variant_interpretation_docker
#        RuntimeAttr? runtime_attr_override
#    }
#
#    Float input_size = size(vcf_file, "GB")
#    Float base_mem_gb = 3.75
#
#    RuntimeAttr default_attr = object {
#                                      mem_gb: base_mem_gb,
#                                      disk_gb: ceil(10 + input_size * 1.5),
#                                      cpu: 1,
#                                      preemptible: 2,
#                                      max_retries: 1,
#                                      boot_disk_gb: 8
#                                  }
#
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#
#    output{
#        File vcf_output = "~{chromosome}.vcf.gz"
#        File no_header_vcf_output = "~{chromosome}.noheader.vcf.gz"
#    }
#
#    command <<<
#        set -euo pipefail
#
#        bcftools index ~{vcf_file}
#
#        bcftools view ~{vcf_file} --regions ~{chromosome} -O z -o  ~{chromosome}.vcf.gz
#
#        bcftools view ~{chromosome}.vcf.gz | grep -v ^## | bgzip -c > ~{chromosome}.noheader.vcf.gz
#
#
#    >>>
#
#    runtime {
#        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
#        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
#        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
#        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
#        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
#        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
#        docker: variant_interpretation_docker
#    }
#}

task mergeGenomicDisorders{
    input{
        Array[File] gd_bed_to_merge
        Array[File] gd_denovo_bed_to_merge
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(gd_bed_to_merge, "GB")
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

    output{
        File gd_output_from_depth = "gd.raw.files.output.ref.txt.gz"
        File gd_denovo_output_from_depth = "gd.denovo.raw.files.output.ref.txt.gz"
    }

    command <<<
        set -euo pipefail

        zcat ~{sep=" " gd_bed_to_merge} | bgzip -c > gd.raw.files.output.ref.txt.gz
       # awk '{print $5"_"$8"_"$9"\t"$6"\t"$7"\t"$8"\t"$9"\t"$1"\t"$2"\t"$3"\t"$4}' gd.raw.files.output.txt | \
       #     bgzip -c > gd.raw.files.output.ref.txt.gz

        zcat ~{sep=" " gd_denovo_bed_to_merge} | bgzip -c > gd.denovo.raw.files.output.ref.txt.gz
       # awk '{print $5"_"$8"_"$9"\t"$6"\t"$7"\t"$8"\t"$9"\t"$1"\t"$2"\t"$3"\t"$4}' gd.raw.files.output.txt | \
       #     bgzip -c > gd.denovo.raw.files.output.ref.txt.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

#task mergeDenovoBedFiles{
#    input{
#        Array[File] bed_files
#        String variant_interpretation_docker
#        RuntimeAttr? runtime_attr_override
#    }
#
#    Float bed_files_size = size(bed_files, "GB")
#    Float base_mem_gb = 3.75
#
#    RuntimeAttr default_attr = object {
#                                      mem_gb: base_mem_gb,
#                                      disk_gb: ceil(10 + (bed_files_size) * 2.0),
#                                      cpu: 1,
#                                      preemptible: 2,
#                                      max_retries: 1,
#                                      boot_disk_gb: 8
#                                  }
#
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#
#    output{
#        File merged_denovo_output = "denovo.merged.bed.gz"
#    }
#
#    command {
#
#        zcat ${bed_files[0]} | head -n+1 > denovo.merged.bed
#        zcat ${sep=" " bed_files} | grep -v ^chrom >> denovo.merged.bed
#        bgzip denovo.merged.bed
#    }
#    runtime {
#        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
#        memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GB"
#        disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
#        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
#        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
#        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
#        docker: variant_interpretation_docker
#    }
#}

task reformatVCF {
    input {
        File vcf
        String docker_path
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible: 0,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File out_bed = "~{prefix}.bed.gz"
        File out_ref_bed = "~{prefix}.ref.bed.gz"
    }

    command <<<
        set -euo pipefail

        echo "Starting svtk"
        svtk vcf2bed -i ALL --include-filters ~{vcf} - | bgzip -c > ~{prefix}.bed.gz
        echo "svtk finished"

        echo "Starting reformat of bed file"
        zcat ~{prefix}.bed.gz | \
            grep -E "DEL|DUP" | \
            awk '{print $1"_"$5"\t"$2"\t"$3"\t"$4"\t"$5}' | \
#            grep -v ^# | \
            bgzip -c > ~{prefix}.ref.bed.gz
        echo "Reformat finished"

  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu, default_attr.cpu])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_path
    preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task getGDvcf {
  input {
    File bed
    File genomic_disorders
    String docker_path
    String prefix
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File out_bed = "~{prefix}.gd.fromVCF.bed.gz"

  }

  command <<<
    set -euo pipefail
    bedtools intersect -wa -wb -f 0.3 -r -a ~{bed} -b ~{genomic_disorders} | bgzip -c > ~{prefix}.gd.fromVCF.bed.gz
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu, default_attr.cpu])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_path
    preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}