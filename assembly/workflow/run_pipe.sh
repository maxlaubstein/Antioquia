#!/bin/bash

#SBATCH -N 1
#SBATCH -t 7-00:00:00
#SBATCH --mem 16G
#SBATCH -n 1
#SBATCH --partition intermediate
#SBATCH -J snakemake

source /n/holylfs05/LABS/informatics/Users/dkhost/mamba/bin/activate snakemake

snakemake --use-conda --rerun-incomplete --profile ../profiles/slurm
