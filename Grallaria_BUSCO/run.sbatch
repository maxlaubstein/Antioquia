#!/bin/bash
#SBATCH --job-name=BUSCO_grallaria
#SBATCH --partition=shared,edwards
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --output=busco_grallaria.%j.log
#SBATCH --error=busco_grallaria.%j.err

source ~/.bashrc 
mamba activate busco_env


busco -i /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Grallaria_repeats/05_full_out/grallaria.p_ctg.fasta.full_mask.soft.fasta \
        -m genome \
        -l passeriformes_odb12 \
        -o BUSCO_output_Grallaria \
        --augustus \
        -c 32

