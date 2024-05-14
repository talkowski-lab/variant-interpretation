version 1.0

import "wes-denovo-helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow exportVDStoVCF {
    input {
        String input_vds
        String hail_docker
        String output_vcf_filename
    }
    
    call helpers.getHailMTSize as getInputVDSSize {
        input:
            mt_uri=input_vds,
            hail_docker=hail_docker
    }

    call exportVDS {
        input:
            input_vds=input_vds,
            output_vcf_filename=output_vcf_filename,
            hail_docker=hail_docker,
            input_size=getInputVDSSize.mt_size
    }

    output {
        File output_vcf = exportVDS.output_vcf
    }
}

task exportVDS {
    input {
        String input_vds
        String output_vcf_filename
        String hail_docker
        Float input_size
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
    cat <<EOF > export_vds.py
    import hail as hl
    import os
    import sys

    input_vds = sys.argv[1]
    output_vcf_filename = sys.argv[2]

    hl.init()

    vds = hl.vds.read_vds(input_vds)
    mt = vds.variant_data
    mt = mt.annotate_entries(GT=hl.vds.lgt_to_gt(mt.LGT, mt.LA))

    rows = mt.entries().select('rsid','gvcf_info').key_by('locus', 'alleles')
    mt = mt.annotate_rows(info=rows[mt.row_key].gvcf_info).drop('gvcf_info')
    hl.export_vcf(mt, output_vcf_filename)
    EOF
    python3 export_vds.py ~{input_vds} ~{output_vcf_filename}
    >>>

    output {
        File output_vcf = output_vcf_filename
    }
}
