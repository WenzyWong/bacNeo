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

Additionally, if you want to run `bacNeo --bacp` to predict bacteria-derived neoantigens, you need to manually install [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), which is not incorporated in conda. Open the link above, click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacNeo --bacp`.

## Tutorial

`bacNeo` has five modules, seperated by five mutually exclusive parameters, i.e., `--download-db`, `--bacc`, `--extract-matrix`, `--bach`, and `--bacp`.

### Commands

- `--download-db` would download all databases and references required in `bacc`, `bach`, and `bach`.

- `--bacc` can extract the number of bacterial reads from genome or transcriptome datasets, and output both raw counts and normalized data (CPM and / or abundance). The command will generate bacterial counts / CPM / abundance of your input taxonomic levels (using 'g' and 's' as an example). The outputs would looks like:

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

- `--extract-matrix` would aid in bacterial read count / CPM-normalization / abundance matrix. Make sure that you have already run `bacNeo --bacc` to produce bacterial information per sample.

- `--bach` can predict HLA alleles for each patient sample from genome datasets. 

- `--bacp` can detect bacterial peptides from proteome datasets and predict HLA-peptide affinities based on results from `bach`. The reference should be provided by users. It is recommended to download reference protome from [UniProt](https://www.uniprot.org/). If no proteome dataset is available, `bacp` could also predict potential neoantigens based on bacteria identified by `bacc`. See usage: `bacp -h`.

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
