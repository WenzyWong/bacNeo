# bacNeo-2.0

## Introduction

The tool **bacNeo-2.0** aim to identify potential intratumour bacteria-derived neoantigens. 

To accomplish the whole antigen-discovery processes, **bacNeo-2.0** requires raw genome, transcriptome, and/or proteome data sets.

The steps processing different omics data are independent. Therefore, you could also use these commands and scripts in `./utils` indenpendently to serve your own specific needs.

To successfully call the commends as well as their sub-functions, please add the directory path of this tool to your environment PATH.

## Commands

- `bacc` can extract bacterial reads from genome or transcriptome datasets. See usage: `bacc -h`.

    ```
    Usage: bacc [ -1 FQ1 ] [ -2 FQ2 ] [-m OMICS] [ -g ZIP ] [ -r REF ] [ -o OUT ] [ -t THREADS ] [ -k KRAKENDB ] [-l TAXONOMY]
        -1 Paired-end clean data (R1) in fastq format.
        -m Type of omics data. 'RNA' for transcriptome, 'WES'/'WGS' for genome.
        -g Whether the fastq file is zipped. 'y' for zipped (.fastq.gz format), 'n' for unzipped (.fastq format).
        -r Reference directory path for hisat2 alignment (if omics data is RNA-seq) or bwa alignment (if omics data is WES/WGS).
        -o Output directory path.
        -t Number of threads (Default threads = 16).
        -k Kraken2 database directory path.
        -l Taxonomy level you would like to generate. Taxonomy levels include: 'd' for Domain, 'p' for Phylum, 'c' for Class, 'o'for Order, 'f' for Family, 'g' for Genus, and 's' for Species. Please ensure the character you input level(s) is(are) included above. If you want to calculate bacterial counts in multiple levels, e.g., in genus and species level, you could input the characters one by one, e.g., -l g -l s.
        If you have multiple sample files and want to run this command repeatedly, it is recommended to make independent directory for each sample.
    ```

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
