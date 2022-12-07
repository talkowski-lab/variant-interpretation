# Base script https://portal.firecloud.org/#methods/Talkowsk-SV/Coverage_plot/10/wdl
version 1.0

import "Structs.wdl"
import "RdTestVisualization.wdl" as rdtest
import "runIgvTrioPlots.wdl" as igv_trio

workflow Module09VisualizeTrio{
    input{
        File? Fasta
        File? Fasta_idx
        File? Fasta_dict

        File varfile
        File pedfile
        String? flags
        String prefix
        File? batch_bincov
        File? sample_batches
        Array[File]? medianfile

        Array[String]? pb_list
        Array[String]? fa_list
        Array[String]? mo_list
        Array[File]? pb_cram_list
        Array[File]? pb_crai_list
        Array[File]? fa_cram_list
        Array[File]? fa_crai_list
        Array[File]? mo_cram_list
        Array[File]? mo_crai_list
        File? nested_repeats
        File? simple_repeats
        File? empty_track

        String sv_base_mini_docker
        String sv_pipeline_rdtest_docker
        String igv_docker

        Boolean run_RD 
        Boolean run_IGV

        RuntimeAttr? runtime_attr_override
        RuntimeAttr? runtime_attr_concatinate
        RuntimeAttr? runtime_attr_rdtest
    }
    

    if(run_RD) {
        Array[File] medianfile_ = select_first([medianfile])
        File batch_bincov_ = select_first([batch_bincov])
        File sample_batches_ = select_first([sample_batches])
        String flags_ = select_first([flags])

        call rdtest.RdTestVisualization as RdTest{
            input:
                prefix = prefix,
                medianfile = medianfile_,
                pedfile = pedfile,
                batch_bincov=batch_bincov_,
                bed = varfile,
                sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
                sample_batches = sample_batches_,
                flags = flags_,
                runtime_attr_rdtest=runtime_attr_rdtest

        }
    }

    if (run_IGV) {   
        Array[String] pb_list_ = select_first([pb_list])
        Array[String] fa_list_ = select_first([fa_list])
        Array[String] mo_list_ = select_first([mo_list])
        Array[File] pb_cram_list_ = select_first([pb_cram_list])
        Array[File] fa_cram_list_ = select_first([fa_cram_list])
        Array[File] mo_cram_list_ = select_first([mo_cram_list])
        Array[File] pb_crai_list_ = select_first([pb_crai_list])
        Array[File] fa_crai_list_ = select_first([fa_crai_list])
        Array[File] mo_crai_list_ = select_first([mo_crai_list])
        File Fasta_ = select_first([Fasta])
        File Fasta_idx_ = select_first([Fasta_idx])
        File Fasta_dict_ = select_first([Fasta_dict])
        File nested_repeats_ = select_first([nested_repeats])
        File simple_repeats_ = select_first([simple_repeats])
        File empty_track_ = select_first([empty_track])

        call igv_trio.IGV_all_samples as igv_plots {
            input:
                pb_list = pb_list_,
                fa_list = fa_list_,
                mo_list = mo_list_,
                pb_cram_list = pb_cram_list_,
                pb_crai_list = pb_crai_list_,
                fa_cram_list = fa_cram_list_,
                fa_crai_list = fa_crai_list_,
                mo_cram_list = mo_cram_list_,
                mo_crai_list = mo_crai_list_,
                nested_repeats = nested_repeats_,
                simple_repeats = simple_repeats_,
                empty_track = empty_track_,
                varfile = varfile,
                Fasta = Fasta_,
                Fasta_dict = Fasta_dict_,
                Fasta_idx = Fasta_idx_,
                prefix = prefix,
                sv_base_mini_docker = sv_base_mini_docker,
                igv_docker = igv_docker,
                runtime_attr_override=runtime_attr_override
        }
    }

    if (run_RD && run_IGV) {
        File igv_plots_tar_gz_pe_ = select_first([igv_plots.tar_gz_pe])
        File RdTest_Plots_ = select_first([RdTest.Plots])

        call concatinate_plots{
            input:
                rd_plots = RdTest_Plots,
                igv_plots = igv_plots_tar_gz_pe,
                prefix = prefix,
                varfile = varfile,
                pedfile = pedfile,
                igv_docker = igv_docker,
                runtime_attr_concatinate = runtime_attr_concatinate
        }
    }

    output{
        File output_plots = select_first([concatinate_plots.plots, RdTest.Plots, igv_plots.tar_gz_pe])
        
    }
}

task concatinate_plots{
    input{
        File rd_plots
        File igv_plots
        String? prefix
        File varfile
        File pedfile
        String igv_docker
        RuntimeAttr? runtime_attr_concatinate
    }

    RuntimeAttr default_attr = object {
        cpu: 1, 
        mem_gb: 7.5,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible: 3,
        max_retries: 1
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
        cat ~{varfile} | cut -f1-5,97 > updated_varfile.bed
        tail -n+2 updated_varfile.bed > ~{varfile}.noheader
        echo 'test'
        python3 /src/MakeRDtest.py \
            ~{varfile}.noheader \
            ~{pedfile} \
            ~{prefix} \
            10000000 \
            ~{prefix}_igv_plots \
            ~{prefix}_rd_plots/ \
            ~{prefix}_igv_rdtest_plots
        tar -czf ~{prefix}_igv_rdtest_plots.tar.gz ~{prefix}_igv_rdtest_plots
    >>>

    output{
        File plots = "~{prefix}_igv_rdtest_plots.tar.gz"
    }


}