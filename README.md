# bacNeo-2.0

## Introduction

The tool **bacNeo** is for identifying intratumor bacteria-derived neoantigens. 

To accomplish the whole antigen-discovery processes, **bacNeo** requires genome, transcriptome, and/or proteome data sets. 

The processes handling different omics data are independent. Therefore, you could also use indenpendent commands in our repository, i.e., `bacc`, `bach` and `bacp`.

To successfully call the commends as well as their sub-functions, the directory path of this tool **must be** added to your environment PATH.

## Commands

- `bacc` can extract bacterial reads from genome or transcriptome datasets. See usage: `bacc -h`.

- `bach` can predict HLA alleles for each patient sample from genome datasets. See usage: `bach -h`.

- `bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. If no proteome dataset available, `bacp` could also use 10-mer winder sliding to chop up reference proteome of interested species, and predict HLA-peptide affinities as well. See usage: `bacp -h`.

## Requirments

### For `bacc`

- If the input data is RNA-seq data: [hisat2](https://daehwankimlab.github.io/hisat2/)

- If the input data is WES/WGS data: [bwa](https://github.com/lh3/bwa)

- [samtools](https://www.htslib.org/)

- [Kraken 2](https://ccb.jhu.edu/data/kraken2_protocol/)

- python (3.7)

### For `bach`

- [bwa](https://github.com/lh3/bwa)

- [samtools](https://www.htslib.org/)

- [hlascan](https://github.com/SyntekabioTools/HLAscan/) (Optional. You could also install hlascan within `bach` by specifying the installation directory using `-s`.)

### For `bacp`

- If proteome datasets are provided: [maxquant](https://anaconda.org/bioconda/maxquant)

- [netMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
