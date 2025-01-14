# bacNeo

## Introduction

The tool **bacNeo** aims to identify potential intratumour bacteria-derived neoantigens. 

To accomplish the whole antigen-discovery processes, **bacNeo** requires raw genome, transcriptome, and/or proteome data sets.

The steps for processing sequencing data are modular and independent. Therefore, you could use `bacc` and `bach` commands, as well as other scripts in `./utils` independently to meet your specific needs.

To successfully call the commends, please add the directory path of this tool to your environment `PATH`.

## Installation

For simplicity, you could clone `bacNeo` from GitHub and run it directory without installing.

First, clone this repository:

```bash
git clone https://github.com/WenzyWong/bacNeo.git
cd bacNeo
```

Then, create an [Anaconda](https://docs.anaconda.com/anaconda/install/) environment with prerequisites using the `bacNeo.yml` file:

```bash
conda env create -n bacNeo -f bacNeo.yml
conda activate bacNeo
```

It is recommended to add the `./bin` path of cloned repository into your `PATH` variable through:

```bash
export PATH=$PATH:/place/of/bacNeo/bin
```

You could test `bacNeo` through running:

```bash
bacNeo -h
```

Additionally, if you want to run `bacp`, you need to manually install [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), which is not incorporated in conda. Click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacp`.

## Tutorial

### Command `bacNeo`

`bacNeo` has two modules, seperated by two mutually exclusive parameters, i.e., `--download-db` and `--extract-matrix`. 

- `--download-db` would download all databases and references required in `bacc`, `bach`, and `bach`.

- `--extract-matrix` would aid in bacterial read count / CPM-normalization / abundance matrix. Make sure that you have already run `bacc` to produce bacterial read counts per sample.

For detailed usage: `bacNeo -h`.

```
Usage: bacNeo [ --download-db DB -t THREADS | --extract-matrix -d DIR -l LEVEL -m METHOD ]
    --download-db       The path you would like to put reference databases into
    -t                  Number of threads (Default threads = 16)
    --extract-matrix    Extract matrix using specified parameters
    -d, --dir           Directory path for matrix extraction
    -l, --level         Taxonomic level for calculation
    -m, --method        Normalization method name, including 'raw_count', 'CPM', and 'abundance'
    -h, --help          Show this help message
```

### Command `bacc`

`bacc` can extract the number of bacterial reads from genome or transcriptome datasets. See usage: `bacc -h`.

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
    Make sure that you have already run bacNeo2 --download-db previously to generate required databases.
    If you have multiple sample files and want to run this command repeatedly, it is recommended to make independent directory for each sample.
```

The command will generate bacterial counts / CPM / abundance of your input taxonomic levels (using 'g' and 's' as an example). The outputs would looks like:

```
bacc/
└── SAMPLE
       ├── counts_g.txt
       ├── counts_s.txt
       ├── normalized_g.txt
       ├── normalized_s.txt
       ├── SAMPLE.bai
       ├── SAMPLE.mpa
       ├── SAMPLE.KRAKEN
       ├── SAMPLE_sorted.bam
       ├── SAMPLE_unmap.bam
       ├── SAMPLE_unmap_R1.fq
       └── SAMPLE_unmap_R2.fq
```

### Command `bach`

`bach` can predict HLA alleles for each patient sample from genome datasets. See usage: `bach -h`.

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
    -d Database directory path
    -g The name(s) of HLA type(s)
    HLA types include: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, MICA, MICB, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5, TAP1, and TAP2
    We recommend you use HLA class I types (A, B, and C), if your are interested in intra-tumour bacterial neoantigens
    If you would like to impute multiple HLA types at once, you could input the types one by one, e.g., -g HLA-A -g HLA-B
    -o Output directory path
    
    Optional:
    -t Number of threads (Default: 16)

    Note: The alignment workflow (-1, -2, -r) and pre-processed BAM workflow (-c) are mutually exclusive
    You must choose one workflow or the other
```

### Command `bacp`

`bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. It is recommended to download reference protome from [UniProt](https://www.uniprot.org/). If no proteome dataset is available, `bacp` could also predict potential neoantigens based on bacteria identified by `bacc`. See usage: `bacp -h`.

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
    1. The proteome workflow (-p) and predicted peptide workflow (non -p) are mutually exclusive
    2. bach must be run first to identify HLA alleles for each patient
    3. For multiple samples, use independent directories for each sample
```

### Example script

Here's an example written in `bash` script, running all the commands above to generate predicted bacterial neoantigens based on WGS data:

```bash
#!/bin/bash
# Replace the following variants
dir_dt="${DIRECTORY_PATH_OF_YOUR_DATA}"
ref="${PATH_OF_YOUR_REFERENCE_FASTA}/hg38.fa"
tool="${DIRECTORY_PATH_YOU_WOULD_LIKE_TO_INSTALL_REQUIRED_TOOLS_INTO}"
db="${DIRECTORY_PATH_YOU_WOULD_LIKE_TO_INSTALL_REQUIRED_DATABASES_INTO}"
out="${DIRECTORY_PATH_OF_OUTPUT}"

# Assuming you have multiple samples
ls "${dir_dt}" | while read sample
do
    # Using paired-end WES fastq data as examples
    fq1="${dir_dt}/${sample}/${sample}_R1.fq.gz"
    fq2="${dir_dt}/${sample}/${sample}_R2.fq.gz"
    
    # Run bacc
    mkdir -p "${out}/bacc"
    out_bacc="${out}/bacc/${sample}"
    mkdir -p "${out_bacc}"
    bacNeo --download-db "${db}"
    bacc -1 "${fq1}" -2 "${fq2}" -m WGS -g -r "${ref}" -o "${out_bacc}" -d "${db}" -l s

    # Run bach, using aligned .bam file from bacc result
    mkdir -p "${out}/bach"
    out_bach="${out}/bach/${sample}"
    mkdir -p "${out_bach}"
    bach -c "${out_bacc}" -s "${tool}" -d "${db}" -g "HLA-A" -o "${out_bach}"

    # Run bacp, using outputs from bacc and bach as inputs
    mkdir -p "${out}/bacp"
    out_bacp="${out}/bacp/${sample}"
    mkdir -p "${out_bacp}"
    bacp -i "${out_bacc}" -d ${db} -a "${out_bach}" -o "${out_bacp}"
done
```

## Requirements

### For `bacNeo`, `bacc` , and `bach`

All requirements are listed in `bacNeo.yml`, and could be installed through the commands listed in the section [Installation](#installation).

### For `bacp`

- `Anaconda` based requirements: install the requirements listed in `bacNeo.yml`, following the the commands in the section [Installation](#installation).

- [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/): You need to manually install it because it is not incorporated in `Anaconda`. Open the link above, click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacp`.
