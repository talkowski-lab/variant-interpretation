# Variant interpretation

This repository contains utilities for consolidating, annotating and filtering structural variants (SVs) from [GATK-SV](https://github.com/broadinstitute/gatk-sv).

## Table of Contents
* [Requirements](#requirements)
* [De novo SV discovery](#denovo)
* [Contact and credits](#contact)

## <a name="requirements">Requirements</a>
### Deployment and execution:
* A [Google Cloud](https://cloud.google.com/) account.
* A workflow execution system supporting the [Workflow Description Language](https://openwdl.org/) (WDL), such as:
  * [Cromwell](https://github.com/broadinstitute/cromwell) (v36 or higher). A dedicated server is highly recommended.
* [cromshell](https://github.com/broadinstitute/cromshell) for interacting with a dedicated Cromwell server.

### Data:
* BED
* VCF
* Pedigree file
* Disorder input
* Raw directory

## <a name="denovo">De novo SV discovery</a>
The main scripts to run this analysis is `runDeNovoSVs.wdl`.    

### Execution
```
> git clone https://github.com/talkowski-lab/variant-interpretation.git
> cd variant-interpretation/wdl
> zip dependencies.zip *
> cromshell submit runDeNovoSVs.wdl /path/to/denovo-svs.json /path/to/config.json dependencies.zip
```

## <a name="contact">Contact and credits</a>
Copyright (c) 2022 Talkowski Lab and The Broad Institute of M.I.T. and Harvard  
Contact: [Alba Sanchis-Juan](mailto:asanchis-juan@mgh.harvard.edu)

Variant interpretation team: Alba Sanchis-Juan, Emma Pierce-Hoffman, Harrison Brand, Mark Walker, Chelsea Lowther, Elise Valkanas