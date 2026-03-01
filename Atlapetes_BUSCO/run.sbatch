#!/bin/bash
#SBATCH --job-name=BUSCO_blancae
#SBATCH --partition=shared,edwards
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --output=busco_blancae.%j.log
#SBATCH --error=busco_blancae.%j.err

source ~/.bashrc 
mamba activate busco_env


busco -i /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Atlapetes_repeats/05_full_out/blancae.p_ctg.fasta.full_mask.soft.fasta \
        -m genome \
        -l passeriformes_odb12 \
        -o BUSCO_output_Atlapetes \
        --augustus \
        -c 32
