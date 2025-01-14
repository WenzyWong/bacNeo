# bacNeo

## Introduction

The tool **bacNeo** aims to identify potential intratumour bacteria-derived neoantigens. 

To accomplish the whole antigen-discovery processes, **bacNeo** requires raw genome, transcriptome, and/or proteome data sets.

The steps for processing sequencing data are modular and independent. Therefore, you could use `bacc` and `bach` commands, as well as other scripts in `./utils` independently to meet your specific needs.

To successfully call the commends, please add the directory path of this tool to your environment `PATH`.

## Installation

For simplicity, you could clone `bacNeo` from GitHub and run it directory without installing.

First, clone this repository:

```
git clone https://github.com/WenzyWong/bacNeo.git
cd bacNeo
```

Then, create an [Anaconda](https://docs.anaconda.com/anaconda/install/) environment with prerequisites using the `bacNeo.yml` file:

```
conda env create -n bacNeo -f bacNeo.yml
conda activate bacNeo
```

It is recommended to add the `./bin` path of cloned repository into your `${PATH}` variable through:

```
export PATH=$PATH:/place/of/bacNeo/bin
```

You could test `bacNeo` through running:

```
bacNeo -h
```

Additionally, if you want to run `bacp`, you need to manually install [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), which is not incorporated in conda. Click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacp`.

## Commands

- `bacNeo` has two modules, seperated by two mutually exclusive parameters, i.e., `--download-db` and `--extract-matrix`. 

    - `--download-db` would download all databases and references required in `bacc`, `bach`, and `bach`.

    - `--extract-matrix` would aid in bacterial read count / CPM-normalization / abundance matrix. Make sure that you have already run `bacc` to produce bacterial read counts per sample.

    For detailed usage: `bacNeo -h`.

    ```
        Usage: bacNeo2 [ --download-db DB | --extract-matrix -d DIR -l LEVEL -m METHOD ]
            --download-db     The path you would like to put reference databases into.
            --extract-matrix  Extract matrix using specified parameters
            -d, --dir        Directory path for matrix extraction
            -l, --level      Taxonomic level for calculation
            -m, --method     Normalization method name
            -h, --help       Show this help message.
    ```

- `bacc` can extract the number of bacterial reads from genome or transcriptome datasets. See usage: `bacc -h`.

    ```
    Usage: bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -g ] [ -r REF ] [ -o OUT ] [ -t THREADS ] [ -d DB ] [-l TAXONOMY] [ -c ]
        -1 Paired-end clean data (R1) in fastq format.
        -2 Paired-end clean data (R2) in fastq format.
        -m Type of omics data. 'RNA' for transcriptome, 'WES'/'WGS' for genome.
        -g [Optional] If the fastq file is zipped.
        -r Reference directory path for hisat2 alignment (if omics data is RNA-seq) or bwa alignment (if omics data is WES/WGS).
        -o Output directory path.
        -t Number of threads (Default threads = 16).
        -d The directory path of your datasebases generated through bacNeo2 --download-db.
        -l Taxonomy level you would like to generate. 
            Taxonomy levels include: 'd' for Domain, 'p' for Phylum, 'c' for Class, 'o'for Order, 'f' for Family, 'g' for Genus, and 's' for Species. 
            Please ensure the level(s) you input is(are) included above. 
            If you would like to calculate bacterial counts in multiple levels, you could input the characters one by one, e.g., -l g -l s.
        -c [Optional] Check the quality and potential contamination of bacterial reads.
        Make sure that you've already run bacNeo2 --download-db previously to generate required databases.
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
    Usage: 
        For alignment workflow:
        bach -1 FQ1 -2 FQ2 -r REF -s SCAN -d DB -g GENES -o OUT [-t THREADS]
        
        For pre-processed BAM workflow (if you have already used WES / WGS mode to run bacc, or you only have .bam file):
        bach -c BACC_PATH -s SCAN -d DB -g GENES -o OUT [-t THREADS]

        Required for alignment workflow:
        -1 Paired-end clean data (R1) in fastq format
        -2 Paired-end clean data (R2) in fastq format
        -r Reference fasta file for bwa alignment, either hg38 or hg19

        Required for pre-processed BAM workflow:
        -c Directory path containing pre-processed BAM files

        Required for both workflows:
        -s The directory path which you would like to install hla-scan into
        -d Database directory path. It's in 'reference' our tool package, named 'HLA-ALL.IMGT'
        -g The name(s) of HLA type(s)
        HLA types include: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, MICA, MICB, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5, TAP1, and TAP2
        We recommend you use HLA class I types (A, B, and C), if your are interested in intra-tumour bacterial neoantigens
        If you would like to impute multiple HLA types at once, you could input the types one by one, e.g., -g HLA-A -g HLA-B
        -o Output directory path
        
        Optional:
        -t Number of threads (Default: 16)

        Note: The alignment workflow (-1, -2, -r) and pre-processed BAM workflow (-c) are mutually exclusive. 
        You must choose one workflow or the other.
    ```

- `bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. It is recommended to download reference protome from [UniProt](https://www.uniprot.org/). If no proteome dataset is available, `bacp` could also predict potential neoantigens based on bacteria identified by `bacc`. See usage: `bacp -h`.

    ```
    Usage: 
        For proteome workflow:
            bacp -p -i MS_DATA -r REF_DIR -d DB_DIR -a ALLELE_DIR -o OUTPUT [-t THREADS]
        
        For predicted peptide workflow:
            bacp -i BACC_DIR -d DB_DIR -a ALLELE_DIR -o OUTPUT [-t THREADS]

        Required for proteome workflow:
        -p             Flag for proteome data analysis
        -i MS_DATA     Directory containing MS data in '.d' format
        -r REF_DIR     Directory containing reference proteome fasta files from UniProt
                        (https://www.uniprot.org/)
                        Recommend using a clean directory with only fasta files

        Required for predicted peptide workflow:
        -i BACC_DIR    Directory containing bacc output

        Required for both workflows:
        -d DB_DIR      Directory containing databases generated by bacNeo2 --download-db
        -a ALLELE_DIR  Directory containing bach results, with sample-specific folders
        -o OUTPUT      Output directory path

        Optional:
        -t THREADS     Number of threads (Default: 16)

        Notes: 
        1. The proteome workflow (-p) and peptide workflow (non-p) are mutually exclusive
        2. Bach must be run first to identify HLA alleles for each patient
        3. For multiple samples, use independent directories for each sample
    ```

## Requirements

### For `bacc` and `bach`

All requirements are listed in `bacNeo.yml`, and could be installed through the commands listed in the section [Installation](#installation).

### For `bacp`

- [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
