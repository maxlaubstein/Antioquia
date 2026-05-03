#!/bin/bash
#SBATCH --job-name=compleasm
#SBATCH --partition=shared,edwards
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH --output=compleasm.%j.log
#SBATCH --error=compleasm.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate compleasm

assembly="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa"
compleasm run  -a $assembly -o compleasm_out  -l aves_odb12 -t 16
