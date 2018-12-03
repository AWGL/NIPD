# Non-Invasive Prenatal Diagnosis Bioinformatics Pipeline

## Introduction

## Installation

The pipeline is designed to run on the Linux (Centos) and Mac OS operating systems.

The pipeline is designed to be run within a Conda environment. Complete the following steps to install all dependencies:

1) Install Miniconda from (https://conda.io/miniconda.html).

2) Clone and enter repository:

```
git clone https://github.com/AWGL/NIPD.git
cd NIPD

```

3) Create a conda environment:

`conda env create -f envs/main.yaml`

4) Activate the environment:

`source activate nipd`

## Configuration

The pipeline is configured using a YAML file. By convention these files are stored in the config/ directory.

Example configuration files are provided.

The input for the pipeline should be in the following structure:

```
Snakefile
config
envs
sample1/
	-sample1_S1_L001_R1_001.fastq.gz
	-sample1_S1_L001_R2_001.fastq.gz
	-sample1_S1_L002_R1_001.fastq.gz
	-sample1_S1_L002_R2_001.fastq.gz
	-sample1_S1_L00n_R1_001.fastq.gz
	-sample1_S1_L00n_R2_001.fastq.gz
sample2/
	-sample2_S2_L001_R1_001.fastq.gz
	-sample2_S2_L001_R2_001.fastq.gz
	-sample2_S2_L002_R1_001.fastq.gz
	-sample2_S2_L002_R2_001.fastq.gz
	-sample2_S2_L00n_R1_001.fastq.gz
	-sample2_S2_L00n_R2_001.fastq.gz
```


## Run

To run locally cd into the directory containing the script and enter the following command:

`snakemake`

To run on cluster enter the following command:

`snakemake --jobs 1 --cluster "qsub -d /share/data/results/snakemake_test/" --directory /share/data/results/snakemake_test/ -s /share/data/results/snakemake_test/Snakefile `

## Output

## References