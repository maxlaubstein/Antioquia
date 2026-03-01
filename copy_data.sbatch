#!/bin/bash
#SBATCH --job-name=Antioquia_data_copy
#SBATCH --partition=sapphire
#SBATCH --time=10:00:00
#SBATCH --mem=4G
#SBATCH --output=copy_pacbio_%j.log
#SBATCH --ntasks=1

#Just copying the raw PacBio reads over from the sequencing facility (Bauer Core_
# Copy PacBio reads
rsync -avh --progress /n/ngsdata/pacbio/251201_REV0053/r84147_20251201_185556/1_B01 ./
