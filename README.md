# COSMIC database preprocessing

This program uses the cosmic VCF file (Coding or Non coding) and Mutant export data (tsv) to generate data of prevalance of variants across tumor types.
### _<span style="color:red">This repo does **NOT** provide access to COSMIC database or license information. Users must obtain their own license to use the COSMIC database.</span>_

#### Software requirements:

- Python v3.8 or higher
- progressbar2 (python package; `pip install progressbar2`)
- sh (python package; `pip install sh`)
- Mac OSX (>=10.11) and Linux (Kernel version 5.x). It has not been tested natively on Windows OS but can be executed inside WSL (Windows Subsystem for Linux)

Hardware requirements:

- 2 CPU cores
- 6-8 GB RAM

Commandline arguments

```
$> python3 process_cosmic_data.py CosmicCodingMuts.vcf.gz CosmicMutantExport.tsv.gz
```
INPUT
- CosmicCodingMut.vcf.gz - VCF file of all coding mutations in the current release. There is a corresponding non-coding VCF file that can be used as input.
- CosmicMutantExport.tsv.gz - A tab separated table of all COSMIC coding point mutations from targeted and genome wide screens from the current release.
- Further details at https://cancer.sanger.ac.uk/cosmic/download

OUTPUT  
Two files are generated
- `<cosmic_vcf_name>_<genome_ver>_processed.txt.gz` - This VCF file can be directly used for generating the Echtvar database. Alternatively, the VCF file can be used for custom annotation pipelines as well.
- `missing_cosmic_ids.txt` - This text file, lists the variants from the COSMIC vcf files that were not found in the `CosmicMutantExport` file.

Example of a variant record in the output VCF file:

```
chr1    685726  .       C       T       32      PASS    COSMIC_ID=COSV60458950;TUMOR_TYPE=thyroid:4|haematopoietic_and_lymphoid_tissue:1
```

##### _NOTE: input files must be gzipped (by default gzipped with downloaded from COSMIC's database)_

