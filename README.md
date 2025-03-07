# bacNeo

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

```
bacNeo path found in: ${THE_PATH_YOU_CLONED_THIS_REPO}
Reference directory: ${THE_PATH_YOU_CLONED_THIS_REPO}/reference


               _                 _   _
              | |               | \ | |
              | |__   __ _  ___ |  \| | ___  ___
              | `_ \ / _` |/ __|| . ` |/ _ \/ _ \
              | |_) | (_| | (__ | |\  |  __/ (_) |
              |____/ \__,_|\___||_| \_|\___|\___/




    Thank you for downloading bacNeo!

    The tool is developed for predicting bacteria-derived neoantigens

    bacNeo contains five modules, please check the usage down below

    Usage:

    --- Module 0 - Download databases
        (Download and build reference databases)
            bacNeo --download-db [-t THREADS]
                --download-db       Download and build reference databases
                -t                  [Optional] Number of threads (Default threads = 16)
        Note: The process is potentially time-consuming

    --- Module 1 - Run bacc
        (Identify bacterial reads and abundance from WGS / WES / RNA-seq data)
            bacNeo --bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -r REF ] [ -o OUT ] [ -t THREADS ] [-l TAXONOMY_LEVEL]

                --bacc              Identify bacterial reads and abundance from WGS / WES / RNA-seq data
                -1                  Paired-end clean data (R1) in fastq format
                -2                  Paired-end clean data (R2) in fastq format
                -m                  Type of omics data. 'RNA' for transcriptome, 'WGS' / 'WES' for genome
                -r                  If '-m RNA' is input: reference directory path for hisat2 alignment
                                    If '-m WGS/WES' is input: reference directory path for bwa alignment
                -o                  Output directory path
                -t                  [Optional] Number of threads (Default threads = 16)
                -l                  Taxonomy level you would like to extract
                                    Taxonomy levels include: 'd' for Domain, 'p' for Phylum, 'c' for Class, 'o'for Order, 'f' for Family, 'g' for Genus, and 's' for Species
                                    Please ensure the level(s) you input is(are) included above
                                    If you would like to calculate bacterial counts and normalized counts in multiple levels, you could input the characters one by one, e.g., -l g -l s
        Notes:
            1. Make sure that you have already run 'bacNeo --download-db' to generate required databases
            2. For multiple samples, it is recommended to use independent directories for each sample

    --- Module 2 - Extract matrix
        (Extract counts / normalized counts after running 'Module 1 - BACC')
            bacNeo --extract-matrix [ -d BACC_OUT_DIR ] [ -l TAXONOMY_LEVEL ] [ -n NORM ]

                --extract-matrix    Extract matrix using specified parameters
                -d, --dir           Directory path for matrix extraction
                -l, --level         Taxonomic level for calculation
                -n, --norm          Normalization method name, including 'raw_count', 'CPM', and 'abundance'
        Note: Make sure that you have already run 'bacNeo --bacc' to extract bacterial reads in all samples

    --- Module 3 - Run bach
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
                -t                  [Optional] Number of threads (Default threads = 16)
        Notes:
            1. Make sure that you have already run 'bacNeo --download-db' to generate required databases
            2. If you have genome data, it is recommended to run the second workflow to save space and time
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
                -t                  [Optional] Number of threads (Default threads = 16).
        Notes:
            1. The proteome workflow (-p) and predicted peptide workflow (non -p) are mutually exclusive
            2. For both (-p) and (non -p) workflow, make sure that you have already run 'bacNeo --download-db' and 'bacNeo --bach'
            3. For (non -p) workflow, make sure that you have already run 'bacNeo --bacc'
            4. For multiple samples, you could input them all at one, based on the directory path of 'bacNeo --bach' and / or 'bacNeo --bacc' output

    -h, --help                      Show this help message
```

Additionally, if you want to run `bacNeo --bacp` to predict bacteria-derived neoantigens, you need to manually install [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), which is not incorporated in conda. Open the link above, click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacNeo --bacp`.

## Tutorial

`bacNeo` has five modules, seperated by five mutually exclusive parameters, i.e., `--download-db`, `--bacc`, `--extract-matrix`, `--bach`, and `--bacp`.

### Commands

- `--download-db` would download all databases and references required in `--bacc`, `--bach`, and `--bach`. The required data would be downloaded and installed in the `reference/` folder. 

    See usage:

    ```
    Usage:

    bacNeo --download-db [-t THREADS]
        --download-db       Download and build reference databases
        -t                  [Optional] Number of threads (Default threads = 16)
    Note: The process is potentially time-consuming
    ```

    After successfully running the command, the folder would look like:

    ```
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
    │   │   │   └── prelim_map.txt
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
    ```

- `--bacc` can extract the number of bacterial reads from genome or transcriptome datasets, and output both raw counts and normalized data (CPM and / or abundance). The command will generate bacterial counts / CPM / abundance of your input taxonomic levels.

    See usage:

    ```
    bacNeo --bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -r REF ] [ -o OUT ] [ -t THREADS ] [-l TAXONOMY_LEVEL]

        --bacc              Identify bacterial reads and abundance from WGS / WES / RNA-seq data
        -1                  Paired-end clean data (R1) in fastq format
        -2                  Paired-end clean data (R2) in fastq format
        -m                  Type of omics data. 'RNA' for transcriptome, 'WGS' / 'WES' for genome
        -r                  If '-m RNA' is input: reference directory path for hisat2 alignment
                            If '-m WGS/WES' is input: reference directory path for bwa alignment
        -o                  Output directory path
        -t                  [Optional] Number of threads (Default threads = 16)
        -l                  Taxonomy level you would like to extract
                            Taxonomy levels include: 'd' for Domain, 'p' for Phylum, 'c' for Class, 'o'for Order, 'f' for Family, 'g' for Genus, and 's' for Species
                            Please ensure the level(s) you input is(are) included above
                            If you would like to calculate bacterial counts and normalized counts in multiple levels, you could input the characters one by one, e.g., -l g -l s
    Notes:
        1. Make sure that you have already run 'bacNeo --download-db' to generate required databases
        2. For multiple samples, it is recommended to use independent directories for each sample
    ```

    The outputs would look like (using 'g' and 's' as an example):

    ```
    bacc/
    └── SAMPLE
        │  ├── counts_g.txt
        │  ├── counts_s.txt
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

- `--extract-matrix` would aid in bacterial read count / CPM-normalization / abundance matrix. Make sure that you have already run `bacNeo --bacc` to produce bacterial information per sample. 

    See usage:

    ```
    bacNeo --extract-matrix [ -d BACC_OUT_DIR ] [ -l TAXONOMY_LEVEL ] [ -n NORM ]

        --extract-matrix    Extract matrix using specified parameters
        -d, --dir           Directory path for matrix extraction
        -l, --level         Taxonomic level for calculation
        -n, --norm          Normalization method name, including 'raw_count', 'CPM', and 'abundance'

    Note: Make sure that you have already run 'bacNeo --bacc' to extract bacterial reads in all samples
    ```

    After successfully running the command, the outputs directory of `bacNeo --bacc` would looks like (assuming that you have already followed the aforementioned steps in `bacNeo --bacc` successfully):

    ```
    bacc/
    ├── matrix_abundance_g.txt
    ├── matrix_abundance_s.txt
    ├── Plot_abundance_g.pdf
    ├── Plot_abundance_s.pdf
    └── SAMPLE
        │  ├── counts_g.txt
        │  ├── counts_s.txt
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

- `--bach` can predict HLA alleles for each patient sample from genome datasets. If you use the sampe genome data to run `bacNeo --bacc` previously, you could skip the alignment process to save time. 

    See usage:

    ```
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
        -t                  [Optional] Number of threads (Default threads = 16)
    Notes:
        1. Make sure that you have already run 'bacNeo --download-db' to generate required databases
        2. If you have genome data, it is recommended to run the second workflow to save space and time
        3. For multiple samples, it is recommended to use independent directories for each sample
    ```

    Taking HLA-A as an example, the output folder would look like:

    ```
    bach/
    └── SAMPLE
           ├── align
           │   ├── SAMPLE_chr6.bam
           │   └── SAMPLE_chr6.bam.bai
           └── HLA-A.txt
    ```

- `--bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. It is recommended to download reference protome from [UniProt](https://www.uniprot.org/). If no proteome dataset is available, `bacp` could also predict potential neoantigens based on bacteria identified by `bacc`. 

    See usage:

    ```
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
        -t                  [Optional] Number of threads (Default threads = 16).
    Notes:
        1. The proteome workflow (-p) and predicted peptide workflow (non -p) are mutually exclusive
        2. For both (-p) and (non -p) workflow, make sure that you have already run 'bacNeo --download-db' and 'bacNeo --bach'
        3. For (non -p) workflow, make sure that you have already run 'bacNeo --bacc'
        4. For multiple samples, you could input them all at one, based on the directory path of 'bacNeo --bach' and / or 'bacNeo --bacc' output
    ```

    You could select peptides accroding to `Strong_binders.csv` and `Weak_binders.csv` with both TAP transportation efficiencty and HLA binding affinity as filtering reference. You could also apply the sequence information in `03_amino_acid_processing/` to other filtering methods of your own interests. The output would look like:

    ```
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

### Example script

Here's an example written in `bash` script, running all the commands above to generate predicted bacterial neoantigens based on WGS data:

```bash
#!/bin/bash
# Replace the following variants
dir_dt="${DIRECTORY_PATH_OF_YOUR_DATA}"
ref="${PATH_OF_YOUR_REFERENCE_FASTA}/hg38.fa"
out="${DIRECTORY_PATH_OF_OUTPUT}"
threads="${THREADS}"

bacNeo --download-db -t "${threads}"
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
    bacNeo --bacc -1 "${fq1}" -2 "${fq2}" -m WGS -r "${ref}" -o "${out_bacc}" -l s -l g -t "${threads}"

    # Run bach, using aligned .bam file from bacc result
    out_bach="${out}/bach/${sample}"
    mkdir -p "${out_bach}"
    bacNeo --bach -c "${out_bacc}" -g "HLA-A" -o "${out_bach}" -t "${threads}"
done
# Run bacp, using outputs from bacc and bach as inputs
bacNeo --bacp -i "${out}/bacc" -a "${out}/bach" -o "${out}/bacp" -t "${threads}"

# Extract bacc normalized matrix, using "abundance" in species-level as an example
bacNeo --extract-matrix -d "${out}/bacc" -l s -n abundance

```

## Requirements

### For modules other than `--bacp`

All requirements are listed in `bacNeo.yml`, and could be installed through the commands listed in the section [Installation](#installation).

### For module `--bacp`

- `Anaconda` based requirements: install the requirements listed in `bacNeo.yml`, following the the commands in the section [Installation](#installation).

- [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/): You need to manually install it because it is not incorporated in `Anaconda`. Open the link above, click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacNeo --bacp`.
