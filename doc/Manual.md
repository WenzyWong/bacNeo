# bacNeo Manual

## Table of Contents

- [Introduction](#introduction)

- [Installation](#installation)

- [Requirements](#requirements)

  - [`bash` Version](#bash-version)

  - [netMHCpan-4.1](#netmhcpan-41)

  - [Harkware Recommendation](#hardware-recommendation)

- [Module Tutorial](#module-tutorial)

  - [Reference Databases](#reference-databases)

  - [BACC](#bacc)

  - [Extract Matrix](#extract-matrix)

  - [BACH](#bach)

  - [BACP](#bacp)

- [Other Utils](#other-utils)

  - [Decontamination](#decontamination)
  
- [Example Script](#example-script)

- [Datasets](#datasets)

  - [Pre-constructed Reference Databases](#pre-constructed-reference-databases)

  - [Test Data](#test-data)

## Introduction

The tool **bacNeo** aims to identify potential intratumour bacteria-derived neoantigens.

To accomplish the whole antigen-discovery processes, **bacNeo** requires raw genome, transcriptome, and / or proteome data sets.

The steps for processing sequencing data are modular and independent. Therefore, you could use `bacNeo --bacc` and `bacNeo --bach` commands, as well as other scripts in `./utils` independently to meet your specific needs.

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

This will output the usage of all modules on your screen:

```plain
bacNeo path found in: ${THE_PATH_YOU_CLONED_THIS_REPO}
Reference directory: ${THE_PATH_YOU_CLONED_THIS_REPO}/reference


               _                 _   _
              | |               | \ | |
              | |__   __ _  ___ |  \| | ___  ___
              | `_ \ / _` |/ __|| . ` |/ _ \/ _ \
              | |_) | (_| | (__ | |\  |  __/ (_) |
              |____/ \__,_|\___||_| \_|\___|\___/




    Thank you for downloading bacNeo!
    
    The tool is developped for predicting bacteria-derived neoantigens, more detailed tutorial is available in doc/Manual.md
    
    bacNeo contains five modules, please check the usage down below

    Usage:

    --- Module 0 - Download databases
        (Download and build reference databases)
        bacNeo --download-db [-t THREADS]
        --download-db       Download and build reference databases
        -t                  [Optional] Number of threads (Default: 16)
        Note: As database-construction is potentially time-consuming, you could also download the pre-constructed databases from Synapse (syn66327848 and syn66514464)

    --- Module 1 - Run bacc
        (Identify bacterial reads and abundance from WGS / WES / RNA-seq data)
        bacNeo --bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -r REF ] [ -o OUT ] [ -t THREADS ] [-l TAXONOMY_LEVEL]

        --bacc              Identify bacterial reads and abundance from WGS / WES / RNA-seq data
        -1                  Paired-end clean data (R1) in fastq format
        -2                  Paired-end clean data (R2) in fastq format
        -m                  Type of omics data. 'RNA' for transcriptome, 'WGS' or 'WES' for genome
        -r                  Reference directory path for hisat2/bwa
                            If '-m RNA': for hisat2 alignment
                            If '-m WGS/WES': for bwa alignment
        -o                  Output directory path
        -t                  [Optional] Number of threads (Default: 16)
        -l                  Taxonomy level(s)
                            Available inputs include: 'd'-Domain, 'p'-Phylum, 'c'-Class, 'o'-Order, 'f'-Family, 'g'-Genus, and 's'-Species
                            If you would like to extract multiple levels, you could input the characters one by one, e.g., -l g -l s
    Notes:
    1. Make sure that you have already run 'bacNeo --download-db' or manually downloaded required databases
    2. For multiple samples, it is recommended to use independent directories for each sample

    --- Module 2 - Extract matrix
        (Extract counts / normalized counts after running 'Module 1 - BACC')
        bacNeo --extract-matrix [ -d BACC_OUT_DIR ] [ -l TAXONOMY_LEVEL ] [ -n NORM ]

        --extract-matrix    Extract matrix using specified parameters
        -d, --dir           Directory path for matrix extraction
        -l, --level         Taxonomic level(s) for calculation
        -n, --norm          Normalization methods, including 'raw_count', 'CPM', and 'abundance'
        
    Notes: Make sure that you have already run BACC to extract bacterial reads in all samples

    --- Module 3 - Run bach
        (Predict HLA-alleles from WGS / WES data)
        bacNeo --bach [ -1 FQ1 ] [ -2 FQ2 ] [ -r REF ] [ -g GENES ] [ -o OUT ] [ -t THREADS ]
        (or if you have already used WGS/WES in BACC, you could also run)
        bacNeo --bach [ -c BACC_PATH ] [ -g GENES ] [ -o OUT ] [ -t THREADS ]

        --bach              Predict HLA-alleles
        -1                  Paired-end clean data (R1) in fastq format
        -2                  Paired-end clean data (R2) in fastq format
        -r                  Reference fasta file for bwa alignment, either hg38 or hg19
        -g                  HLA type(s), including: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, MICA, MICB, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5, TAP1, and TAP2
                            If you would like to impute multiple HLA types at once, input them one by one, e.g., -g HLA-A -g HLA-B
        -o                  Output directory path
        -c                  Directory path containing pre-processed BAM files
        -t                  [Optional] Number of threads (Default: 16)
    Notes:
    1. Make sure that you have already run 'bacNeo --download-db' or manually downloaded required databases
    2. If you input genome data, it is recommended to run the second workflow to save space and time
    3. For multiple samples, it is recommended to use independent directories for each sample

    --- Module 4 - Run bacp
        (Predict bacterial neoantigens based on proteome data)
        bacNeo --bacp [ -p ] [ -i MS_INPUT ] [ -r REF_DIR ] [ -a ALLELE_DIR ] [ -o OUTPUT ] [ -t THREADS ]
        (or predict bacterial neoantigens based on previously identified bacterial reads from 'Module 1 - Run bacc')
        bacNeo --bacp [ -i BACC_OUT_DIR ] [ -a ALLELE_DIR ] [ -o OUTPUT ] [ -t THREADS ]

        --bacp              Predict bacterial neoantigens
        -p                  Flag for proteome data analysis, input this parameter only when you have proteome dataset
        -i                  Directory containing MS data in '.d' format (-p)
                            or
                            Directory containing bacc output (non -p)
        -r                  [Only required when '-p' is input]
                            Directory containing reference proteome fasta files from UniProt (https://www.uniprot.org/)
                            Recommend a clean directory with only fasta files
        -a                  Directory containing bach results, with sample-specific folders
        -o                  Output directory path
        -t                  [Optional] Number of threads (Default: 16).
    Notes:
    1. The proteome workflow (-p) and predicted peptide workflow (non -p) are mutually exclusive
    2. For both (-p) and (non -p) workflow, make sure that you have already constructed reference databases and run 'bacNeo --bach'
    3. For (non -p) workflow, make sure that you have already run 'bacNeo --bacc'
    4. For multiple samples, you could input them all at once, based on the directory path of 'bacNeo --bach' and / or 'bacNeo --bacc' output

    -h, --help                      Show this help message
```

Additionally, if you want to run `bacNeo --bacp` to predict bacteria-derived neoantigens, you need to manually install [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), which is not incorporated in conda (see [Requirenments](#requirements)-[netMHCpan-4.1](#netmhcpan-41)).

## Requirements

### `bash` Version

The majority of the codes were written in `bash shell`. The developed version of `bash` was 5.1.16. Make sure the `bash` version in your machine is newer than version 5.0. You can use `echo $BASH_VERSION` to check the version.

If the default `/bin/bash` doesn't meet our requirement and updating the `bash` version may harm your environment, you could also replace `#!/bin/bash` to `#!/${PATH_OF_YOUR_BACNEO_ANACONDA_ENVS}/bin/bash` in `./bin/bacNeo`.

### `tcsh` Version

Running BACP is dependent on netMHCpan, which is written in `tcsh`. The developed version of `tcsh` was 6.18.01. Make sure the `tcsh` version in your machine is newer than version 6.

### netMHCpan-4.1

As `bacNeo --bacp` uses [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), a lisenced software which requires manually application, you need to manually install it. Open the link above, click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacNeo --bacp`.

### Hardware Recommendation

- The pre-contructed reference databases require more than 100 GB of disk spaces. If you download and construct databases via commands, that would cause approximately 350 GB of disk spaces.

- Multi-threading is recommended, and the default thread number is 16. You could also choose running `bacNeo` within one thread or more threads, depending on your hardware circumstances.

- Running `bacNeo` using 16 threads would cause approximately 8 GB of memory, please ensure memory limits of your machine is above 8 GB if you use the default parameters.

## Module Tutorial

`bacNeo` has five modules, three of which are core functional modules. Running `bacNeo` is also seperated by five mutually exclusive parameters, i.e., `--download-db`, `--bacc`, `--extract-matrix`, `--bach`, and `--bacp`.

### Reference Databases

`--download-db` would download all databases and references required in `--bacc`, `--bach`, and `--bach`. The required data would be downloaded and installed in the `reference/` folder. You can also download the pre-contructed databases manually to save time (see [Datasets](#datasets)-[[Pre-constructed Reference Databases](#pre-constructed-reference-databases)]).

See usage:

```plain
Usage:

bacNeo --download-db [-t THREADS]
--download-db       Download and build reference databases
-t                  [Optional] Number of threads (Default: 16)
Note: The process is potentially time-consuming
```

After successfully running the command, the folder would look like:

```plain
reference/
├── bac_na
│   ├── hash.k2d
│   ├── library
│   │   ├── archaea
│   │   │   ├── assembly_summary.txt
│   │   │   ├── library.fna
│   │   │   ├── library.fna.masked
│   │   │   ├── manifest.txt
│   │   │   └── prelim_map.txt
│   │   ├── bacteria
│   │   │   ├── assembly_summary.txt
│   │   │   ├── library.fna
│   │   │   ├── library.fna.masked
│   │   │   ├── manifest.txt
│   │   │   ├── prelim_map.txt
│   │   │   ├── xaa
│   │   │   ├── xab
│   │   │   └── xac
│   │   ├── human
│   │   │   ├── assembly_summary.txt
│   │   │   ├── library.fna
│   │   │   ├── manifest.txt
│   │   │   └── prelim_map.txt
│   │   ├── plasmid
│   │   │   ├── library.fna
│   │   │   ├── library.fna.masked
│   │   │   ├── manifest.txt
│   │   │   └── prelim_map.txt
│   │   ├── UniVec_Core
│   │   │   ├── library.fna
│   │   │   ├── library.fna.masked
│   │   │   ├── prelim_map.txt
│   │   │   └── UniVec_Core
│   │   └── viral
│   │       ├── assembly_summary.txt
│   │       ├── library.fna
│   │       ├── library.fna.masked
│   │       ├── manifest.txt
│   │       └── prelim_map.txt
│   ├── opts.k2d
│   ├── seqid2taxid.map
│   ├── taxo.k2d
│   └── taxonomy
│       ├── accmap.dlflag
│       ├── citations.dmp
│       ├── delnodes.dmp
│       ├── division.dmp
│       ├── gc.prt
│       ├── gencode.dmp
│       ├── images.dmp
│       ├── merged.dmp
│       ├── names.dmp
│       ├── nodes.dmp
│       ├── nucl_gb.accession2taxid
│       ├── nucl_wgs.accession2taxid
│       ├── prelim_map.txt
│       ├── readme.txt
│       ├── taxdump.dlflag
│       ├── taxdump.tar.gz
│       └── taxdump.untarflag
├── CheckM2_database
│   └── uniref100.KO.1.dmnd
├── HLA-ALL.IMGT
└── hla_scan_r_v2.1.4

11 directories, 54 files
```

In `reference/bac_na/`, you could remove all directories and files other than `.k2d` files to save spaces. That would make the structure of `reference/` look like the following one, which is identical to manually downloaded reference (see [Datasets](#datasets) - [Pre-constructed Reference Databases](#pre-constructed-reference-databases)):

```plain
reference/
├── bac_na
│   ├── hash.k2d
│   ├── opts.k2d
│   └── taxo.k2d
├── CheckM2_database
│   └── uniref100.KO.1.dmnd
├── HLA-ALL.IMGT
└── hla_scan_r_v2.1.4

2 directories, 6 files
```

## BACC

`--bacc` can extract the number of bacterial componants from genome or transcriptome datasets, and output raw counts and normalized data (CPM and/or abundance) with charts and graphics to assist downstream analysis. The command will generate sample-specific bacterial counts/CPM/abundance within your input taxonomic level(s).

See usage:

```plain
bacNeo --bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -r REF ] [ -o OUT ] [ -t THREADS ] [-l TAXONOMY_LEVEL]

    --bacc              Identify bacterial reads and abundance from WGS/WES/RNA-seq data
    -1                  Paired-end clean data (R1) in fastq format
    -2                  Paired-end clean data (R2) in fastq format
    -m                  Type of omics data. Use 'RNA' for transcriptome, 'WGS'/'WES' for genome
    -r                  For transcriptome: reference directory path for hisat2 alignment, e.g., ${PATH_OF_HISAT}/hg38
                        For genome: reference directory path for bwa alignment, e.g., ${PATH_OF_BWA}/hg38.fa
    -o                  Output directory path
    -t                  [Optional] Number of threads (Default: 16)
    -l                  Taxonomy level(s) you would like to extract
                        Taxonomy levels include: 'd'-Domain, 'p'-Phylum, 'c'-Class, 'o'-Order, 'f'-Family, 'g'-Genus, and 's'-Species
                        Please ensure the level(s) you input is(are) included above
                        If you would like to calculate bacterial counts and normalized counts in multiple levels, you could input the characters one by one, e.g., -l g -l s
Notes:
1. Make sure that you have already run 'bacNeo --download-db' to generate required databases, or manually download reference databases from Synapse
2. For multiple samples, it is recommended to use independent directories for each sample
```

The outputs would look like (using 'g' and 's' as an example):

```plain
bacc/
└── SAMPLE
    ├── counts_g.txt
    ├── counts_s.txt
    ├── normalized_g.txt
    ├── normalized_s.txt
    ├── SAMPLE.KRAKEN
    ├── SAMPLE.mpa
    ├── SAMPLE_sorted.bam
    ├── SAMPLE_sorted.bam.bai
    ├── SAMPLE.standard
    ├── SAMPLE_unmap.bam
    ├── SAMPLE_unmap_R1.fq
    └── SAMPLE_unmap_R2.fq

2 directories, 12 files
```

### Extract Matrix

`--extract-matrix` would aid in bacterial read count / CPM-normalization / abundance matrix, and generate phlygenetic tree based on species abundant across samples. Make sure that you have already run `bacNeo --bacc` to produce bacterial information per sample.

See usage:

```plain
bacNeo --extract-matrix [ -d BACC_OUT_DIR ] [ -l TAXONOMY_LEVEL ] [ -n NORM ]

    --extract-matrix    Extract matrix using specified parameters
    -d, --dir           Directory path for matrix extraction
    -l, --level         Taxonomic level for calculation
    -n, --norm          Normalization method name, including 'raw_count', 'CPM', and 'abundance'

Note: Make sure that you have already run 'bacNeo --bacc' to extract bacterial reads in all samples
```

After successfully running the command, the outputs directory of `bacNeo --bacc` would looks like (assuming that you have already followed the aforementioned steps in `bacNeo --bacc` successfully):

```plain
bacc/
├── matrix_abundance_g.txt
├── matrix_abundance_s.txt
├── out_500.tree
├── out.otu.csv
├── out.tree
├── out_tree_30.png
├── out_tree_500.nwk
├── out_tree_500.png
├── out_tree.nwk
├── out_tree.png
├── Plot_abundance_g.pdf
├── Plot_abundance_s.pdf
├── standard_reports
│   └── SAMPLE.standard
└── SAMPLE
    │  ├── counts_g.txt
    │  └── counts_s.txt
    ├── normalized_g.txt
    └── normalized_s.txt
       ├── SAMPLE.bai
       ├── SAMPLE.mpa
       ├── SAMPLE.KRAKEN
       ├── SAMPLE_sorted.bam
       ├── SAMPLE_unmap.bam
       ├── SAMPLE_unmap_R1.fq
       └── SAMPLE_unmap_R2.fq
```

### BACH

`--bach` can predict HLA alleles for each patient sample from genome datasets. If you use the sampe genome data to run `bacNeo --bacc` previously, you could skip the alignment process to save time.

See usage:

```plain
(Predict HLA-alleles from WGS / WES data)
    bacNeo --bach [ -1 FQ1 ] [ -2 FQ2 ] [ -r REF ] [ -g GENES ] [ -o OUT ] [ -t THREADS ]
(or if you have already used WGS / WES in 'Module 1 - Run bacc', you could also run)
bacNeo --bach [ -c BACC_PATH ] [ -g GENES ] [ -o OUT ] [ -t THREADS ]

    --bach              Predict HLA-alleles
    -1                  Paired-end clean data (R1) in fastq format
    -2                  Paired-end clean data (R2) in fastq format
    -r                  Reference fasta file for bwa alignment, either hg38 or hg19
    -g                  The name(s) of HLA type(s)
                        HLA types include: HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, MICA, MICB, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5, TAP1, and TAP2
                        It is recommended to use HLA class I types (A, B, and C), if your are interested in intra-tumour bacterial neoantigens
                        If you would like to impute multiple HLA types at once, you could input the types one by one, e.g., -g HLA-A -g HLA-B
    -o                  Output directory path
    -c                  Directory path containing pre-processed BAM files
    -t                  [Optional] Number of threads (Default: 16)
Notes:
    1. Make sure that you have already run 'bacNeo --download-db' to generate required databases
    2. If you have genome data, it is recommended to run the second workflow to save space and time
    3. For multiple samples, it is recommended to use independent directories for each sample
```

Taking HLA-A as an example, the output folder would look like:

```plain
bach/

└── SAMPLE
    ├── align
    │   ├── SAMPLE_chr6.bam
    │   └── SAMPLE_chr6.bam.bai
    └── HLA-A.txt

2 directories, 3 files
```

### BACP

- `--bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. It is recommended to download reference protome from [UniProt](https://www.uniprot.org/). If no proteome dataset is available, `bacp` could also predict potential neoantigens based on bacteria identified by `bacc`.

    See usage:

    ```plain
    (Predict bacterial neoantigens based on proteome data)
        bacNeo --bacp [ -p ] [ -i MS_INPUT ] [ -r REF_DIR ] [ -a ALLELE_DIR ] [ -o OUTPUT ] [ -t THREADS ]
    (or predict bacterial neoantigens based on previously identified bacterial reads from 'Module 1 - Run bacc')
        bacNeo --bacp [ -i BACC_OUT_DIR ] [ -a ALLELE_DIR ] [ -o OUTPUT ] [ -t THREADS ]

        --bacp              Predict bacterial neoantigens
        -p                  Flag for proteome data analysis, input this parameter only when you have proteome dataset
        -i                  Directory containing MS data in '.d' format (-p)
                            or
                            Directory containing bacc output (non -p)
        -r                  [Only required when '-p' is input]
                            Directory containing reference proteome fasta files from UniProt (https://www.uniprot.org/)
                            Recommend a clean directory with only fasta files
        -a                  Directory containing bach results, with sample-specific folders
        -o                  Output directory path
        -t                  [Optional] Number of threads (Default: 16).
    Notes:
        1. The proteome workflow (-p) and predicted peptide workflow (non -p) are mutually exclusive
        2. For both (-p) and (non -p) workflow, make sure that you have already run 'bacNeo --download-db' and 'bacNeo --bach'
        3. For (non -p) workflow, make sure that you have already run 'bacNeo --bacc'
        4. For multiple samples, you could input them all at one, based on the directory path of 'bacNeo --bach' and / or 'bacNeo --bacc' output
    ```

    You could select peptides accroding to `Strong_binders.csv` and `Weak_binders.csv` with both TAP transportation efficiencty and HLA binding affinity as filtering reference. You could also apply the sequence information in `03_amino_acid_processing/` to other filtering methods of your own interests. The output would look like:

    ```plain
    bacp/
    ├── 00_allele_summary
    │   ├── Bar_distribution_alleles.pdf
    │   ├── SAMPLE1.txt
    │   ├── SAMPLE2.txt
    │   └── SAMPLE3.txt
    ├── 01_na_fastas
    │   ├── SAMPLE1_bacreads.fasta
    │   ├── SAMPLE2_bacreads.fasta
    │   └── SAMPLE3_bacreads.fasta
    ├── 02_predicted_protein
    │   ├── SAMPLE1
    │   │   ├── checkm2.log
    │   │   ├── diamond_output
    │   │   │   └── DIAMOND_RESULTS.tsv
    │   │   ├── protein_files
    │   │   │   └── SAMPLE1_bacreads.faa
    │   │   └── quality_report.tsv
    │   ├── SAMPLE2
    │   │   ├── checkm2.log
    │   │   ├── diamond_output
    │   │   │   └── DIAMOND_RESULTS.tsv
    │   │   └── protein_files
    │   │       └── SAMPLE2_bacreads.faa
    │   └── SAMPLE3
    │       ├── checkm2.log
    │       ├── diamond_output
    │       │   └── DIAMOND_RESULTS.tsv
    │       ├── protein_files
    │       │   └── SAMPLE3_bacreads.faa
    │       └── quality_report.tsv
    ├── 03_amino_acid_processing
    │   ├── SAMPLE1_bacpep_cdhit.fa
    │   ├── SAMPLE1_bacpep_cdhit.fa.clstr
    │   ├── SAMPLE1_bacpep_windowed.txt
    │   ├── SAMPLE1_bacpep_windowed.txt.temp
    │   ├── SAMPLE2_bacpep_cdhit.fa
    │   ├── SAMPLE2_bacpep_cdhit.fa.clstr
    │   ├── SAMPLE2_bacpep_windowed.txt
    │   ├── SAMPLE2_bacpep_windowed.txt.temp
    │   ├── SAMPLE3_bacpep_cdhit.fa
    │   ├── SAMPLE3_bacpep_cdhit.fa.clstr
    │   ├── SAMPLE3_bacpep_windowed.txt
    │   └── SAMPLE3_bacpep_windowed.txt.temp
    ├── 04_affinity_with_HLAs
    │   ├── SAMPLE1
    │   │   ├── SAMPLE1_HLA-A02_07.xls
    │   │   └── SAMPLE1_HLA-A33_03.xls
    │   ├── SAMPLE2
    │   │   ├── SAMPLE2_HLA-A02_07.xls
    │   │   └── SAMPLE2_HLA-A11_01.xls
    │   └── SAMPLE3
    │       ├── SAMPLE3_HLA-A02_06.xls
    │       └── SAMPLE3_HLA-A11_01.xls
    ├── Strong_binders.csv
    └── Weak_binders.csv
    ```

## Other utils

### Decontamination

To facilitate more strict detection of tumour microbiome, we also provide user-callable scripts for decontamination. `/utils/create_tree.py` can output files for decontamination while creating phylogenic trees. `/utils/decontam_analysis.R` can then utilise these files to run decontamination processes. Both normal samples and negative control samples are required to run `decontam_analysis.R`. The usage is as follows:

```bash
Usage: Rscript decontam_analysis.R <samples_prefix> <negatives_prefix> [output_prefix]
```

## Example script

Here's an example written in `bash` script, running all the commands above to generate predicted bacterial neoantigens based on RNA-seq clean data:

```bash
#!/bin/bash
# Replace the following variants
dir_dt="${DIRECTORY_PATH_OF_YOUR_DATA}"
ref="${PATH_OF_YOUR_REFERENCE_FASTA}"
out="${DIRECTORY_PATH_OF_OUTPUT}"
threads="${THREADS}"

# Remove the following hash tags if you want to run our test data
#prefetch SRR33242434 --output-directory ${dir_dt}
#prefetch SRR33242436 --output-directory ${dir_dt}
#prefetch SRR33242438 --output-directory ${dir_dt}

#ls "${dt_dir}" | while read id
#do
#    fastq-dump --split-files "${dt_dir}/${id}/${id}.sra" --gzip --outdir "${dt_dir}/${id}"
#    rm "${dt_dir}/${id}/${id}.sra"
#done

# If pre-constructed reference databases have not been downloaded, remove the hash below
# bacNeo --download-db -t "${threads}"
mkdir -p "${out}/bacc"
mkdir -p "${out}/bach"
mkdir -p "${out}/bacp"

ls "${dir_dt}" | while read sample
do
    fq1="${dir_dt}/${sample}/${sample}.R1.fq.gz"
    fq2="${dir_dt}/${sample}/${sample}.R2.fq.gz"
    
    # Run bacc
    out_bacc="${out}/bacc/${sample}"
    mkdir -p "${out_bacc}"
    # Using both genus and species levels as an example
    bacNeo --bacc -1 "${fq1}" -2 "${fq2}" -m RNA -r "${ref}" -o "${out_bacc}" -l s -t "${threads}"

    # Run bach, using aligned .bam file from bacc result
    out_bach="${out}/bach/${sample}"
    mkdir -p "${out_bach}"
    bacNeo --bach -c "${out_bacc}" -g "HLA-A" -o "${out_bach}" -t "${threads}"
done

# Extract bacc normalized matrix, using "abundance" in species-level as an example
bacNeo --extract-matrix -d "${out}/bacc" -l s -n abundance

# Run bacp, using outputs from bacc and bach as inputs
bacNeo --bacp -i "${out}/bacc" -a "${out}/bach" -o "${out}/bacp" -t "${threads}"
```

## Datasets

### Pre-constructed Reference Databases

If you would like to download the pre-constructed reference databases, pleas put them into the `reference/` directory in your `bacNeo` path, and ignore the `bacNeo --download-db` module. Databases are available on Synapse (accession: [syn66327848](https://www.synapse.org/Synapse:syn66327848/files/) and [syn66514464](https://www.synapse.org/Synapse:syn66514464/files/)).

### Test Data

We used the raw WES data from BioProject (accession: [PRJNA1253793](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1253793) - SRR33242434, SRR33242436, and SRR33242438).
