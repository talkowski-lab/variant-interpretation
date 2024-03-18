# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs2.wdl"
import "RdVisualization.wdl" as rdtest
import "CreateIgvCramPlots.wdl" as igv_cram
import "CreateIgvEvidencePlots.wdl" as igv_evidence

workflow VisualizePlots{
    input{
        File varfile
        File pedfile
        String prefix
        File? batch_bincov
        File? sample_batches
        File? batch_medianfile
        File? fam_ids
        Int? max_size


        Int ?igv_max_window
        File? rd_outliers
        File? sample_pe_sr
        File? sample_crai_cram
        String? buffer
        File? reference
        File? reference_index
        Boolean? file_localization
        Boolean? requester_pays
        Boolean? is_snv_indel
#        File? regeno_file

        String sv_base_mini_docker
        String sv_pipeline_rdtest_docker
        String igv_docker
        String variant_interpretation_docker

        Boolean run_RD 
        Boolean run_IGV
        Boolean run_evidence_plots
        Boolean run_cram_plots

        RuntimeAttr? runtime_attr_run_igv
        RuntimeAttr? runtime_attr_igv
        RuntimeAttr? runtime_attr_cpx
        RuntimeAttr? runtime_attr_concatinate
        RuntimeAttr? runtime_attr_rdtest
        RuntimeAttr? runtime_attr_reformat_pe
        RuntimeAttr? runtime_attr_reformat_sr
        RuntimeAttr? runtime_attr_update_pe_sr
    }
    String max_size_ = select_first([max_size, 10000000])

    Boolean is_snv_indel_ = if defined(is_snv_indel) then select_first([is_snv_indel]) else false
#    Boolean is_snv_indel_ = select_first([is_snv_indel])

    #Update complex bed file
    if (is_snv_indel_ == false){
        call updateCpxBed{
            input:
                varfile = varfile,
                variant_interpretation_docker = variant_interpretation_docker,
                runtime_attr_override = runtime_attr_cpx
        }
    }

    #creates RD plots for DELs and DUPs
    if(run_RD) {
        File batch_medianfile_ = select_first([batch_medianfile])
        File batch_bincov_ = select_first([batch_bincov])
        File sample_batches_ = select_first([sample_batches])
        File rd_outliers_ = select_first([rd_outliers])
#        File regeno_file_ = select_first([regeno_file])

        call rdtest.RdTestVisualization as RdTest{
            input:
                prefix = prefix,
                ped_file = pedfile,
                fam_ids = fam_ids,
                batch_medianfile = batch_medianfile_,
                batch_bincov=batch_bincov_,
                bed = select_first([updateCpxBed.bed_output, varfile]),
#                regeno=regeno_file_,
                sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                variant_interpretation_docker = variant_interpretation_docker,
                outlier_samples = rd_outliers_,
                sample_batches = sample_batches_,
                runtime_attr_rdtest=runtime_attr_rdtest

        }
    }

    #creates IGV plots for all variants (proband will be the top plot if it has affected status = 2 in ped file)
    if (run_IGV) {   
        File buffer_ = select_first([buffer,500])
        File reference_ = select_first([reference])
        File reference_index_ = select_first([reference_index])
        Boolean file_localization_ = if defined(file_localization) then select_first([file_localization]) else false
        if(run_evidence_plots){
            File sample_pe_sr_ = select_first([sample_pe_sr])
            call igv_evidence.IGV_all_samples as igv_evidence_plots {
                input:
                    ped_file = pedfile,
                    sample_pe_sr = sample_pe_sr_,
                    buffer = buffer_,
                    fam_ids = fam_ids,
                    varfile = select_first([updateCpxBed.bed_output, varfile]),
                    reference = reference_,
                    reference_index = reference_index_,
                    prefix = prefix,
                    is_snv_indel = is_snv_indel_,
                    file_localization = file_localization_,
                    sv_base_mini_docker = sv_base_mini_docker,
                    igv_docker = igv_docker,
                    variant_interpretation_docker = variant_interpretation_docker,
                    runtime_attr_run_igv = runtime_attr_run_igv,
                    runtime_attr_igv = runtime_attr_igv,
                    runtime_attr_cpx = runtime_attr_cpx,
                    runtime_attr_reformat_pe = runtime_attr_reformat_pe,
                    runtime_attr_reformat_sr = runtime_attr_reformat_sr,
                    runtime_attr_update_pe_sr = runtime_attr_update_pe_sr   
            }
        }
        
        if(run_cram_plots){
            File sample_crai_cram_ = select_first([sample_crai_cram])
            Int igv_max_window_ = if defined(igv_max_window) then select_first([igv_max_window]) else 150000
            Boolean requester_pays_ = if defined(requester_pays) then select_first([requester_pays]) else false
            call igv_cram.IGV_all_samples as igv_cram_plots {
                input:
                    ped_file = pedfile,
                    sample_crai_cram = sample_crai_cram_,
                    buffer = buffer_,
                    fam_ids = fam_ids,
                    varfile = select_first([updateCpxBed.bed_output, varfile]),
                    igv_max_window = igv_max_window_,
                    reference = reference_,
                    file_localization = file_localization_,
                    requester_pays = requester_pays_,
                    is_snv_indel = is_snv_indel_,
                    reference_index = reference_index_,
                    prefix = prefix,
                    sv_base_mini_docker = sv_base_mini_docker,
                    igv_docker = igv_docker,
                    variant_interpretation_docker = variant_interpretation_docker,
                    runtime_attr_run_igv = runtime_attr_run_igv,
                    runtime_attr_igv = runtime_attr_igv,
                    runtime_attr_cpx = runtime_attr_cpx
            }
        }
    }

    #creates a concatinated image with the IGV plot as the top pane and the RD plot as the bottom pane
    if (run_RD && run_IGV) {
        File igv_plots_tar_gz_pe_ = select_first([igv_cram_plots.tar_gz_pe, igv_evidence_plots.tar_gz_pe])
        File RdTest_Plots_ = select_first([RdTest.Plots])

        call concatinate_plots{
            input:
                rd_plots = RdTest_Plots_,
                igv_plots = igv_plots_tar_gz_pe_,
                prefix = prefix,
                varfile = select_first([updateCpxBed.bed_output, varfile]),
                pedfile = pedfile,
                igv_docker = igv_docker,
                runtime_attr_concatinate = runtime_attr_concatinate
        }
    }

    output{
        File output_plots = select_first([concatinate_plots.plots, RdTest.Plots, igv_evidence_plots.tar_gz_pe, igv_cram_plots.tar_gz_pe])
        
    }
}

task concatinate_plots{
    ##To do: take out as taks in now in ConcatinatePlots.wdl
    input{
        File rd_plots
        File igv_plots
        String? prefix
        File varfile
        File pedfile
        String igv_docker
        Int max_size_
        RuntimeAttr? runtime_attr_concatinate
    }

    Float input_size = size(select_all([rd_plots, igv_plots, varfile, pedfile]), "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: 10,
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_concatinate, default_attr])

    

    runtime {
        cpu: select_first([runtime_attr.cpu, default_attr.cpu])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: igv_docker
        preemptible: select_first([runtime_attr.preemptible, default_attr.preemptible])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    command <<<
        set -eu -o pipefail    

        tar -zxf ~{rd_plots}
        tar -zxf ~{igv_plots}
        mkdir ~{prefix}_igv_rdtest_plots
        cat ~{varfile} | gunzip | cut -f1-6 > updated_varfile.bed
        tail -n+2 updated_varfile.bed > ~{varfile}.noheader
        echo 'test'
        python3 /src/MakeRDtest.py \
            ~{varfile}.noheader \
            ~{pedfile} \
            ~{prefix} \
            ~{max_size_} \
            ~{prefix}_igv_plots \
            ~{prefix}_rd_plots/ \
            ~{prefix}_igv_rdtest_plots
        tar -czf ~{prefix}_igv_rdtest_plots.tar.gz ~{prefix}_igv_rdtest_plots
    >>>

    output{
        File plots = "~{prefix}_igv_rdtest_plots.tar.gz"
    }


}

task updateCpxBed{
    input{
        File varfile
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(varfile, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
                                      mem_gb: base_mem_gb,
                                      disk_gb: ceil(10 + input_size * 1.5),
                                      cpu: 1,
                                      preemptible: 2,
                                      max_retries: 1,
                                      boot_disk_gb: 8
                                  }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bed_output = "final.denovo.merged.cpx_split.bed.gz"
    }

    command <<<
        set -euo pipefail

        Rscript /src/variant-interpretation/scripts/updateCpx.R ~{varfile}
        bgzip final.denovo.merged.cpx_split.bed
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