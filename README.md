# Variant interpretation

This repository contains utilities for consolidating, annotating and filtering structural variants (SVs) from [GATK-SV](https://github.com/broadinstitute/gatk-sv).

## Table of Contents
* [Requirements](#requirements)
* [De novo SV discovery](#denovo)
* [Contact and credits](#contact)

## <a name="requirements">Requirements</a>
### Deployment and execution:
* A [Google Cloud](https://cloud.google.com/) account.
* A workflow execution system supporting the [Workflow Description Language](https://openwdl.org/) (WDL), such as [Cromwell](https://github.com/broadinstitute/cromwell) (v36 or higher). A dedicated server is highly recommended.
* [cromshell](https://github.com/broadinstitute/cromshell) for interacting with a dedicated Cromwell server.
 
## <a name="denovo">De novo SV discovery</a>
Filters <i>de novo</i> SV obtained from [GATK-SV](https://github.com/broadinstitute/gatk-sv) to improve specificity. The filtering is based on a combination of site-specific and sample-specific metrics:
* Frequency filter
  * gnomAD AF <= 0.01 OR in genomic disorder (GD)
  * Other parents in the cohort: AF <= 0.01 OR AC <= 3 
* Absent call in both parents
* By type specific filters:
  * CPX,CTX,INV: keep all
  * INS: require raw evidence support
  * DEL,DUP: require raw evidence support AND,
    * Large (>5KB) CNVs require depth support AND no parent overlap (50%) AND no RDCN evidence in parent
    * Small (<1Kb) CNVs with SR-only evidence require input on both sides

### Input data:
* VCF file: annotated with [GATK-SV AnnotateVcf](https://github.com/broadinstitute/gatk-sv#annotatevcf-in-development).
* BED file: obtained from the main VCF GATK-SV file and converted with `svtk vcf2bed --i all --include-filters vcf bed`. 
* PED file: with columns `family_id, subject_id, paternal_id, maternal_id, sex, affected, family_size`.
* Disorder calls: list of SV names that overlap with genomic disorder regions.
* Directory with raw calls from [GATK-SV GatherBatchEvidence](https://github.com/broadinstitute/gatk-sv#gatherbatchevidence) module with the following columns: `chrom, start, end, svtype, sample`
  * The `chrom` column corresponds to a string formed by `chrom_svtype_sample` so the overlap can be performed at a sample and SV type level.

### Execution
The main scripts to run this analysis is `runDeNovoSVs.wdl`.
```
> git clone https://github.com/talkowski-lab/variant-interpretation.git
> cd variant-interpretation/wdl
> zip dependencies.zip *
> cromshell submit runDeNovoSVs.wdl /path/to/denovo-svs.json /path/to/config.json dependencies.zip
```

## <a name="contact">Contact and credits</a>
Copyright (c) 2022 Talkowski Lab and The Broad Institute of M.I.T. and Harvard  
Contact: [Alba Sanchis-Juan](mailto:asanchis-juan@mgh.harvard.edu)

Variant interpretation team: Alba Sanchis-Juan, Harrison Brand, Emma Pierce-Hoffman, Mark Walker, Chelsea Lowther, Elise Valkanas.