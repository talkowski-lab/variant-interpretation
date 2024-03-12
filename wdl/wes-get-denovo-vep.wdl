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

workflow getDenovoVEPandLOEUF {
    input {
        File loeuf_file
        String filtered_mt
        String denovo_ht
        String hail_docker
    }

    call helpers.getHailMTSize as getHailMTSize {
        input:
        mt_uri=filtered_mt,
        hail_docker=hail_docker
    }

    call getDenovoVEP {
        input:
        filtered_mt=filtered_mt,
        denovo_ht=denovo_ht,
        hail_docker=hail_docker,
        input_size=getHailMTSize.mt_size
    }

    call annotateLOEUF {
        input:
        de_novo_vep=getDenovoVEP.de_novo_vep,
        loeuf_file=loeuf_file,
        hail_docker=hail_docker
    }

    output {
        File de_novo_vep = annotateLOEUF.de_novo_vep_loeuf
    }
}

task getDenovoVEP {
    input {
        String filtered_mt
        String denovo_ht
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

    String prefix = basename(filtered_mt, '_wes_denovo_basic_filtering.mt')
    command <<<
    cat <<EOF > get_denovo_vep.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    filtered_mt = sys.argv[1]
    denovo_ht = sys.argv[2]
    cores = sys.argv[3]
    mem = int(np.floor(float(sys.argv[4])))

    prefix = os.path.basename(filtered_mt).split('_wes_denovo_basic_filtering.mt')[0]

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{mem}g",
                    "spark.driver.cores": cores,
                    "spark.driver.memory": f"{mem}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.read_matrix_table(filtered_mt)
    de_novo_results = hl.read_table(denovo_ht)
    mt = mt.semi_join_rows(de_novo_results.key_by('locus', 'alleles'))
    df = mt.rows().to_pandas()
    df.to_csv(f"{prefix}_wes_final_denovo_vep.txt", sep='\t', index=False)
    EOF

    python3 get_denovo_vep.py ~{filtered_mt} ~{denovo_ht} ~{cpu_cores} ~{memory}
    >>>

    output {
        File de_novo_vep = "~{prefix}_wes_final_denovo_vep.txt"
    }
}

task annotateLOEUF {
    input {
        File de_novo_vep
        File loeuf_file
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size([de_novo_vep, loeuf_file], 'GB')
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

    String prefix = basename(de_novo_vep, '_wes_final_denovo_vep.txt')
    command <<<
    cat <<EOF > annotate_loeuf.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os
    import ast

    de_novo_vep = sys.argv[1]
    loeuf_file = sys.argv[2]

    prefix = os.path.basename(de_novo_vep).split('_wes_final_denovo_vep.txt')[0]

    vep_res = pd.read_csv(de_novo_vep, sep='\t')
    vep_res['alleles'] = vep_res.alleles.apply(ast.literal_eval)
    vep_res['ID'] = vep_res.locus + ':' + vep_res.alleles.str.join(':')
    vep_res.columns = vep_res.columns.str.replace('info.', '')

    loeuf = pd.read_csv(loeuf_file, sep='\t')
    loeuf.index = loeuf.gene_name

    def get_genes_csq(csq):
        genes = []
        for ind_csq in csq:
            gene = ind_csq.split('|')[3]
            if gene!='':
                genes.append(gene)
        return list(set(genes))

    vep_res['CSQ'] = vep_res.CSQ.replace({np.nan: "[]"}).apply(ast.literal_eval)
    vep_res['all_genes'] = vep_res.CSQ.apply(get_genes_csq)

    all_genes = vep_res.all_genes.apply(pd.Series).stack().unique()

    loeuf_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF'].to_dict()
    loeuf_tile_vals = loeuf.loc[np.intersect1d(loeuf.index, all_genes), 'LOEUF_tile'].to_dict()

    vep_res['LOEUF'] = vep_res.all_genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_vals).min())
    vep_res['LOEUF_tile'] = vep_res.all_genes.apply(lambda gene_list: pd.Series(gene_list).map(loeuf_tile_vals).min())

    vep_res.to_csv(f"{prefix}_wes_final_denovo_vep_loeuf.txt", sep='\t', index=False)
    EOF

    python3 annotate_loeuf.py ~{de_novo_vep} ~{loeuf_file}
    >>>

    output {
        File de_novo_vep_loeuf = "~{prefix}_wes_final_denovo_vep_loeuf.txt"
    }
}