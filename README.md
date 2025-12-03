# Integrating Gene Expression, Mutation and Copy Number Data to Identify Driver Genes of Recurrent Chromosome-Arm Losses

This repository contains the code associated with the paper **“Integrating Gene Expression, Mutation and Copy Number Data to Identify Driver Genes of Recurrent Chromosome-Arm Losses.”** The manuscript is available [here](https://doi.org/10.1016/j.celrep.2025.116455).

---

## Introduction

The code in this repository is designed to identify driver genes responsible for cancer type–specific recurring arm losses using data from 20 cancer types provided by TCGA. It is written in R (tested on R versions 4.3 and 4.4) and uses [Snakemake](https://snakemake.readthedocs.io/) to orchestrate most of the execution pipeline. Additionally, **MutSig2CV** and **GISTIC2** are employed.

---

## Setup

- The required R packages and Linux dependencies are listed in the `dependencies.txt` file.
- Instructions for installing Snakemake are available [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Both GISTIC2 and MutSig2CV can be installed directly or run via Docker containers (see below).

---

## Running the Code

The execution of the code involves several steps:

1. **Running the Snakemake Pipeline**  
   Calculates the mutation and CNV rates given on arm loss/no arm loss for each cancer type and frequently lost arm. It also performs differential gene expression and pathway analysis and prepares the `.seg` and mutation files for the following step.

2. **Running MutSig2CV and GISTIC2**  
   These tools are executed for each pair of cancer type and frequently lost arm on the 2 groups of samples separately  (with and without arm loss).

3. **Running PRODIGY on GISTIC2 Results**  
   This step (also using Snakemake) applies the PRODIGY algorithm to the GISTIC2 results.

4. **Summarizing the Results**

5. **Plotting (Optional)**

---

### Running the Snakemake Pipeline

The configuration file that specifies the data, logs, and output directories is located at `snakemake_scripts/config.yaml`. By default, these directories are created under the root of the repository. The pipeline expects the following inputs:

- **GISTIC2 Results:**  
  Cancer type–specific results downloaded from [Broad GDAC](https://gdac.broadinstitute.org/).

- **Raw Gene Counts:**  
  Downloaded using TCGAbiolinks (htseq-count data from the legacy archive, which is no longer available).

- **FPKM Gene Counts:**  
  Downloaded using TCGAbiolinks.

Since the TCGA legacy archive is no longer accessible through TCGAbiolinks, and to simplify data access, the required data has been uploaded to [TODO](url). Please place this data in the `data` directory.

To run the pipeline, execute the following commands:

```bash
cd snakemake_scripts
snakemake all --cores NUMBER_OF_CORES_TO_USE
```

---

### Running MutSig2CV and GISTIC2

We encountered errors when trying to run multiple instances of these tools in parallel, so we created several Docker containers for each tool. The Docker image for **GISTIC2** is available [here](https://hub.docker.com/r/shixiangwang/gistic) and for **MutSigCV** [here](https://hub.docker.com/r/genepattern/mutsigcv).

Note that the MutSigCV container does not actually have MutSig2CV installed and is only used for its Matlab setup. Alternatively, you can use a Docker image that includes MutSig2CV.

Download MutSig2CV by executing the following commands from the root of the repository:

```bash
wget http://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSig2CV.tar.gz
tar -xvzf MutSig2CV.tar.gz
```

The markers, reference, and CNV files required for running **GISTIC2** are also available [here](TODO) and should be placed in the `data` directory. These files were also downloaded from [Broad GDAC](https://gdac.broadinstitute.org/).

There are 460 runs of each tool, so the process can take a long time. Using 5 containers for each tool, the runs took approximately 72 hours. They can be executed in parallel. Note that for step 3, only the GISTIC2 results are needed.

To run all GISTIC instances, first insert the container names in the `container_names` list in `code/run_gistic.sh`, then run the script from the `code` directory:

```bash
cd code
./run_gistic.sh out data log/GISTIC
```

To run all MutSig2CV instances, first insert the container names in the `container_names` list in `code/run_mutsig.sh`, then run the script from the `code` directory (preferably in a different terminal to run these in parallel):

```bash
cd code
./run_mutsig.sh out data log/MutSig
```
