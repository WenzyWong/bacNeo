# bacNeo

## Introduction

![BioDraft](doc/BioDraft.png)

The tool **bacNeo** aims to identify potential intratumour bacteria-derived neoantigens.

To accomplish the whole antigen-discovery processes, **bacNeo** requires raw genome, transcriptome, and / or proteome data sets.

The steps for processing sequencing data are modular and independent. Therefore, you could use `bacNeo --bacc` and `bacNeo --bach` commands, as well as other scripts in `./utils` independently to meet your specific needs.

To successfully call the commends, please add the directory path of this tool to your environment `PATH`.

## User manual

Please refer to user manual in [`doc/Manual.md`](doc/Manual.md) for detailed guidance.

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

## Requirements

### `bash` Version

The majority of the codes were written in `bash shell`. The developed version of `bash` was 5.1.16. Make sure the `bash` version in your machine is newer than version 5.0. You can use `echo $BASH_VERSION` to check the version.

If the default `/bin/bash` doesn't meet our requirement and updating the `bash` version may harm your environment, you could also replace `#!/bin/bash` to `#!/${PATH_OF_YOUR_BACNEO_ANACONDA_ENVS}/bin/bash` in `./bin/bacNeo`.

### `tcsh` Version

Running BACP is dependent on netMHCpan, which is written in `tcsh`. The developed version of `tcsh` was 6.18.01. Make sure the `tcsh` version in your machine is newer than version 6.

### netMHCpan-4.1

As `bacNeo --bacp` uses [netMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/), a lisenced software which requires manually application, you need to manually install it. Open the link above, click `Downloads`, and choose `Version 4.1b` - `Linux`. After filling in and submitting the form, you could download and install it to successfully run `bacNeo --bacp`.

## Datasets

### Pre-constructed reference databases

If you would like to download the pre-constructed reference databases, pleas put them into the `reference/` directory in your `bacNeo` path, and ignore the `bacNeo --download-db` module. Databases are available on Synapse ([syn66327848](https://www.synapse.org/Synapse:syn66327848/files/) and [syn66514464](https://www.synapse.org/Synapse:syn66514464/files/)).

### Test data

We used three pre-treatment RNA-seq datasets randomized from BioProject (accession: [PRJEB23709](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP105482&o=bytes_l%3Ad) - ERR2208968, ERR2208969, and ERR2208946).
