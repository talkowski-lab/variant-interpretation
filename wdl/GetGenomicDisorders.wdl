version 1.0
    
import "Structs.wdl"
import "ReformatRawFiles.wdl" as raw
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "DeNovoSVsScatter.wdl" as runDeNovo
import "SplitVcf.wdl" as getBatchedVcf

workflow deNovoSV {

    input {
        File ped_input
        File python_config
        File vcf_file
        Array[String] contigs
        File genomic_disorder_input
        File batch_raw_file
        File batch_depth_raw_file
        File exclude_regions
        File sample_batches
        File batch_bincov_index
        Int records_per_shard
        String prefix
        File? fam_ids
        String variant_interpretation_docker
        String sv_pipeline_updates_docker
        RuntimeAttr? runtime_attr_gd
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_raw_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_merge_bed
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_raw_divide_by_chrom
        RuntimeAttr? runtime_attr_raw_reformat_bed
        RuntimeAttr? runtime_attr_merge_final_bed_files
        RuntimeAttr? runtime_attr_create_plots
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_attr_clean_ped
        RuntimeAttr? runtime_attr_call_outliers
        RuntimeAttr? runtime_attr_get_batched_files
        RuntimeAttr? runtime_attr_merge_gd
        RuntimeAttr? runtime_attr_batch_vcf
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_attr_merge
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

        call getBatchedVcf.getBatchedVcf as getBatchedVcf {
            input:
                vcf_file = vcf_file,
                prefix = prefix,
                samples = getBatchedFiles.samples,
                records_per_shard = 10000,
                variant_interpretation_docker = variant_interpretation_docker,
                sv_pipeline_updates_docker = sv_pipeline_updates_docker,
                runtime_attr_batch_vcf = runtime_attr_batch_vcf,
                runtime_override_shard_vcf = runtime_override_shard_vcf,
                runtime_attr_merge = runtime_attr_merge
        }

    }

    #makes a ped file of singletons, duos, and trios for input into the de novo script (only including families of interest)
    call cleanPed{
        input:
            ped_input = ped_input,
            vcf_input = select_first([getBatchedVcf.split_vcf, vcf_file]),
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_clean_ped
    }

   
    #splits raw files into probands and parents and reformats to have chrom_svtype_sample as the first column for probands and chrom_svtype_famid as the first column for parents
    call raw.reformatRawFiles as reformatDepthRawFiles {
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
        #generates a list of genomic disorder regions in the vcf input as well as in the depth raw files
        call getGenomicDisorders{
            input:
                genomic_disorder_input=genomic_disorder_input,
                ped = cleanPed.cleaned_ped,
                vcf_file = select_first([getBatchedVcf.split_vcf, vcf_file]),
                depth_raw_file_proband = reformatDepthRawFiles.reformatted_proband_raw_files[i],
                depth_raw_file_parents = reformatDepthRawFiles.reformatted_parents_raw_files[i],
                chromosome=contigs[i],
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_gd
        }

    
    
    #merges the genomic disorder region output from each chromosome to compile a list of genomic disorder regions
    call mergeGenomicDisorders{
        input:
            genomic_disorder_input=getGenomicDisorders.gd_output_from_depth_raw_files,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_gd
    }

    output {
        File cleaned_ped = cleanPed.cleaned_ped
        File final_denovo_nonOutliers = callOutliers.final_denovo_nonOutliers_output
        File final_denovo_outliers = callOutliers.final_denovo_outliers_output
        File final_denovo_nonOutliers_plots = createPlots.output_plots
        Array [File] denovo_output_annotated = getDeNovo.per_chromosome_annotation_output_file
        File gd_depth = mergeGenomicDisorders.gd_output_from_depth
        File gd_vcf = getGenomicDisorders.gd_output_from_final_vcf[1]
        
    }
}


task getGenomicDisorders{
    input{
        File vcf_file
        File ped
        File depth_raw_file_proband
        File depth_raw_file_parents
        String chromosome
        File genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([vcf_file, ped, genomic_disorder_input, depth_raw_file_parents, depth_raw_file_proband]), "GB")
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
        File gd_output_from_final_vcf = "gd.variants.from.final.vcf.txt.gz"
        File gd_output_from_depth_raw_files = "~{chromosome}.gd.variants.in.depth.raw.files.txt.gz"
        File gd_output_for_denovo = "annotated.gd.variants.names.txt"
    }

    command <<<
        set -euxo pipefail

        bedtools intersect -wa -wb -f 0.3 -r -a ~{vcf_file} -b ~{genomic_disorder_input} | cut -f 3 |sort -u > annotated.gd.variants.names.txt
        
        echo "Done with first line"

        bedtools intersect -wa -wb -f 0.3 -r -a ~{vcf_file} -b ~{genomic_disorder_input} > gd.variants.from.final.vcf.txt
        bgzip gd.variants.from.final.vcf.txt

        echo "Done with GD from vcf"
        
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
        bedtools intersect -wa -wb -a ~{chromosome}.coverage.proband.list.txt -b sorted.depth.proband.bed.gz | sort | uniq | awk -v OFS="\t" '{print $1,$2,$3,$4,$9,$10,$11,$12,$13}' > ~{chromosome}.coverage.proband.txt

        # IMPORTANT "bedtools coverage" throws an error when sorted.depth.proband.bed.gz contains a sample, but the query file "gd.per.sample.txt" (so the ped file) doesn't contain the corresponding sample. 
        echo "done with coverage in depth variants". This error fix should be done in prior steps.

        cat ~{chromosome}.coverage.parents.txt ~{chromosome}.coverage.proband.txt > ~{chromosome}.coverage.txt

        echo "done with cat"
        bedtools intersect -wa -wb -f 0.3 -sorted -a gd.per.sample.txt -b sorted.depth.proband.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.proband.no.r.txt
        bedtools intersect -wa -wb -f 0.3 -sorted -a gd.per.family.txt -b sorted.depth.parents.bed.gz > ~{chromosome}.gd.variants.in.depth.raw.file.parents.no.r.txt
      
        echo "done with intersect no -r"

        cat ~{chromosome}.gd.variants.in.depth.raw.file.proband.no.r.txt ~{chromosome}.gd.variants.in.depth.raw.file.parents.no.r.txt > ~{chromosome}.remove.txt

        echo "done with cat"

        bedtools intersect -v -wb -b ~{chromosome}.remove.txt -a ~{chromosome}.coverage.txt > ~{chromosome}.kept.coverage.txt

        echo "done with grep"
        # to remove sample/family tag from chr in final bed file   
        cat ~{chromosome}.gd.variants.in.depth.raw.file.proband.txt \
            ~{chromosome}.gd.variants.in.depth.raw.file.parents.txt \
            ~{chromosome}.kept.coverage.txt |\
            awk -v OFS="\t" '{sub(/_.*/, "", $1); print}' |\
            awk -v OFS="\t" '{sub(/_.*/, "", $5); print}' | sort | uniq > ~{chromosome}.gd.variants.in.depth.raw.files.txt
        bgzip ~{chromosome}.gd.variants.in.depth.raw.files.txt
        echo "done with cat"
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

task mergeGenomicDisorders{
    input{
        Array[File] genomic_disorder_input
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(genomic_disorder_input, "GB")
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
        File gd_output_from_depth = "gd.raw.files.output.txt.gz"
    }

    command {
        set -euo pipefail

        zcat ${sep=" " genomic_disorder_input} > gd.raw.files.output.txt
        bgzip gd.raw.files.output.txt
    }

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
        File cleaned_ped = "subset_cleaned_ped.txt"
    }

    command {

        Rscript /src/variant-interpretation/scripts/cleanPed.R ${ped_input}
        cut -f2 cleaned_ped.txt | tail -n+2 > all_samples.txt
        bcftools query -l ${vcf_input} > samples_to_include_in_ped.txt
        grep -w -v -f samples_to_include_in_ped.txt all_samples.txt > excluded_samples.txt
        grep -w -f excluded_samples.txt cleaned_ped.txt | cut -f1 | sort -u > excluded_families.txt
        grep -w -v -f excluded_families.txt cleaned_ped.txt > subset_cleaned_ped.txt

    }

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
