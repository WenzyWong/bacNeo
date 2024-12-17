# bacNeo-2.0

## Introduction

The tool **bacNeo-2.0** aims to identify potential intratumour bacteria-derived neoantigens. 

To accomplish the whole antigen-discovery processes, **bacNeo-2.0** requires raw genome, transcriptome, and/or proteome data sets.

The steps for processing sequencing data are modular and independent. Therefore, you could use `bacc` and `bach` commands, as well as other scripts in `./utils` independently to meet your specific needs.

To successfully call the commends, please add the directory path of this tool to your environment `PATH`.

## Installation

For simplicity, you could clone `bacNeo-2.0` from GitHub and run it directory without installing.

First, clone this repository:

```
git clone https://github.com/WenzyWong/bacNeo-2.0.git
cd bacNeo-2.0
```

Then, create an [Anaconda](https://docs.anaconda.com/anaconda/install/) environment with prerequisites using the `bacNeo.yml` file:

```
conda env create -n bacNeo -f bacNeo.yml
conda activate bacNeo
```

It is recommended to add the path of cloned repository into your `${PATH}` variable through:

```
export PATH=$PATH:/place/of/bacNeo-2.0
```

You could test `bacNeo-2.0` through running:

```
bacc -h
```

Additionally, if you want to run `bacp`, you need to manually install [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), which is not incorporated in conda. Click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacp`.

## Commands

- `bacc` can extract the number of bacterial reads from genome or transcriptome datasets. See usage: `bacc -h`.

    ```
    Usage: bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -g ZIP ] [ -r REF ] [ -o OUT ] [ -t THREADS ] [ -k KRAKENDB ] [-l TAXONOMY]
        -1 Paired-end clean data (R1) in fastq format.
        -2 Paired-end clean data (R2) in fastq format.
        -m Type of omics data. 'RNA' for transcriptome, 'WES'/'WGS' for genome.
        -g Whether the fastq file is zipped. 'y' for zipped (.fastq.gz format), 'n' for unzipped (.fastq format).
        -r Reference directory path for hisat2 alignment (if omics data is RNA-seq) or bwa alignment (if omics data is WES / WGS).
        -o Output directory path.
        -t Number of threads (Default threads = 16).
        -k Kraken2 database directory path.
        -l Taxonomy level you would like to generate. 
           Taxonomy levels include: 'd' for Domain, 'p' for Phylum, 'c' for Class, 'o'for Order, 'f' for Family, 'g' for Genus, and 's' for Species. 
           Please ensure the level(s) you input is(are) included above. 
           If you would like to calculate bacterial counts in multiple levels, you could input the characters one by one, e.g., -l g -l s.
        If you have multiple sample files and want to run this command repeatedly, it is recommended to make independent directory for each sample.
    ```

    Test `bacc` using the RNA-seq data in `./testdata/RNA-seq` with the following command (replace `${HISAT}` and `${DB}` with your own path for hisat2 reference and the place you would like to create a bacterial reference database):

    ```
    bacc -1 testdata/RNA-seq/T001_R1.fq.gz -2 testdata/RNA-seq/T001_R2.fq.gz -m RNA -g y -r ${HISAT} -o output/bacc -k ${DB} -l g -l s
    ```
    
    The command will generate bacterial counts in genus and species levels in `./output/bacc`. The outputs should be identical to:

    ```
    bacc/
    └── T001
        ├── counts_g.txt
        ├── counts_s.txt
        ├── T001.bai
        ├── T001.mpa
        ├── T001.o
        ├── T001_sorted.bam
        ├── T001_unmap.bam
        ├── T001_unmap_R1.fq
        └── T001_unmap_R2.fq
    ```
    
- `bach` can predict HLA alleles for each patient sample from genome datasets. See usage: `bach -h`.

    ```
    Usage: bach [ -1 FQ1 ] [ -2 FQ2 ] [ -r REF ] [ -s SCAN ] [ -d DB ] [ -g GENES ] [ -o OUT ] [ -t THREADS ]
        -1 Paired-end clean data (R1) in fastq format.
        -2 Paired-end clean data (R2) in fastq format.
        -r Reference fasta file for bwa alignment, either hg38 or hg19.
        -s The directory path which you would like to install hla-scan into.
        -d Database directory path. It's in 'reference' our tool package, named 'HLA-ALL.IMGT'.
        -g The name(s) of HLA type(s).
            HLA types include: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, MICA, MICB, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5, TAP1, and TAP2.
            We recommend you use HLA class I types (A, B, and C), if your are interested in intra-tumour bacterial neoantigens.
            If you would like to impute multiple HLA types at once, you could input the types one by one, e.g., -g HLA-A -g HLA-B.
        -o Output directory path.
        -t Number of threads (Default threads = 16).
        If you have multiple sample files and want to run this command repeatedly, it is recommended to make independent directory for each sample.
    ```

- `bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. If no proteome dataset available, `bacp` could also use high abundant species from the results of `bacc`, and predict HLA-peptide affinities as well. See usage: `bacp -h`.

## Requirements

### For `bacc` and `bach`

All requirements are listed in `bacNeo.yml`, and could be installed through the commands listed in the section [Installation](#installation).

### For `bacp`

- [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
