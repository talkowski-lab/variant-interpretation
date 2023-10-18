# Variant interpretation

This repository contains utilities for consolidating, annotating and filtering structural variants (SVs) from [GATK-SV](https://github.com/broadinstitute/gatk-sv).

## Table of Contents
* [Requirements](#requirements)
* [Family filtering](#familyfiltering)
* [De novo SV discovery](#denovo)
* [Contact and credits](#contact)

## <a name="requirements">Requirements</a>
### Deployment and execution:
* A [Google Cloud](https://cloud.google.com/) account.
* A workflow execution system supporting the [Workflow Description Language](https://openwdl.org/) (WDL), such as [Cromwell](https://github.com/broadinstitute/cromwell) (v36 or higher). A dedicated server is highly recommended.
* [cromshell](https://github.com/broadinstitute/cromshell) for interacting with a dedicated Cromwell server or a Terra account
* A vcf file output from the [GATK-SV pipeline](https://github.com/broadinstitute/gatk-sv)
 
## <a name="familyfiltering">Annotation and Filtering of Possible Clinically Relevant SVs</a>
Additional annotations and filtering are performed for each family given a input pedigree file. The input is a joint-genotyped VCF obtained from [GATK-SV](https://github.com/broadinstitute/gatk-sv) and the output is a reformated BED file.

### Filtering
* AF <= 0.3 and gnomAD_V2_AF <= 0.01 OR SV overlaps a genomic disorder region (explained below).

### Annotation
* Genomic disorders: all SVs overlapping a list of genomic disorders regions are flagged in the final output. Currenly it requires a minimum reciprocarl overlap with the GD of 30%.
* Constraint: extra columns are added for SVs involving genes that are constrained for LoF or missense variantsion (using pLI, pRec and EO metrics).
* Gene list: a colum that determines if each SV involves a gene present in the gene list provided is included (Note: this is currently using HGNC Symbols).
* OMIM: a colum that determines if each SV involves a gene present in OMIM is included (Note: this is currently using HGNC Symbols).
* Mode of Inheritance:
    * Absent in unaffected: a variant is flagged if it's present in unaffected individuals, and present only in affected. This is particularly relevant for trios and other family structures to identify <i>de novo</i> variants or inherited from affected parents (Column name: `FILT_ABSENT_UNAFF`)
    * Multiple SVs in same gene: a variant is flagged if there are multiple variants in the same gene for the same individual. If trio is available, it only returns if they are in trans (Column name: `MULT_HIT`)
    * Autosomal recessive: a variant is flagged if genotype for affected individual/s is `1/1` and is not for unaffected (Column name: `FILT_AR`)
    * X-linked: a variant is flagged as `FILT_XLR` column if they are in male affected individual/s and they are `1/1` or `1` genotypes.
    * Inherited: possible intereting inherited variants from unaffected parents are flagged in the `FILT_INHERITED` column if they are in the input gene list AND the cohort allele count <= 10 AND the number of unaffected indivuals in the cohort is <= 5.

## <a name="denovo">De novo SV discovery</a>
Filters <i>de novo</i> SV obtained from [GATK-SV](https://github.com/broadinstitute/gatk-sv) to improve specificity of <i>de novo</i> SV calls. The filtering is based on a combination of site-specific and sample-specific metrics:
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
* PED file: with columns `FamID, IndividualID, FatherID, MotherID, Gender, Affected`
* Disorder input: bed file of SV sites that are genomic disorder regions
* Directory with raw calls from [GATK-SV GatherBatchEvidence](https://github.com/broadinstitute/gatk-sv#gatherbatchevidence) module with the following columns: `chrom, start, end, svtype`
* Raw files: txt file with first column as the batches and second column as the raw files generated from [GATK-SV ClusterBatch] (https://github.com/broadinstitute/gatk-sv/blob/main/wdl/ClusterBatch.wdl). This list must contain the raw files from all callers except for depth
* Depth raw files: txt file with first column as the batches and second column as the depth raw files generated from [GATK-SV ClusterBatch](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/ClusterBatch.wdl). This list must contain the raw files from only the depth caller
* Blacklist regions: bed file with a list of regions you want to exclude from the de novo analysis with the following columns: `chrom, start, end`
* Coverage matrices: txt file with first column as batch name, the second column as coverage matrix generated from [GATK-SV GatherBatchEvidence](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GatherBatchEvidence.wdl), and third column as coverage matrix index file generated from [GATK-SV GatherBatchEvidence](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GatherBatchEvidence.wdl)
* Sample and batch input: txt file with samples in first column and their repsective batches in second column
* Python config file: a python file defining input parameters

### Execution
The main scripts to run this analysis is `runDeNovoSVs.wdl`, which can be run on Terra or with Cromshell with the following commands:
```
> git clone https://github.com/talkowski-lab/variant-interpretation.git
> cd variant-interpretation/wdl
> zip dependencies.zip *
> cromshell submit runDeNovoSVs.wdl /path/to/denovo-svs.json /path/to/config.json dependencies.zip
```

## <a name="contact">Contact and credits</a>
Copyright (c) 2022 Talkowski Lab and The Broad Institute of M.I.T. and Harvard  
Contact: [Alba Sanchis-Juan](mailto:asanchis-juan@mgh.harvard.edu) and [Nicole Calamari](mailto:ncalamari@mgh.harvard.edu)

Variant interpretation team: Alba Sanchis-Juan, Harrison Brand, Emma Pierce-Hoffman, Mark Walker, Chelsea Lowther, Elise Valkanas.
