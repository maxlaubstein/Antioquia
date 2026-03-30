#!/bin/bash
#SBATCH --job-name=scaffolding
#SBATCH --partition=shared,edwards
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --output=scaffolding.%j.log
#SBATCH --error=scaffolding.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16


source ~/.bashrc
mamba activate chromap_env

assembly="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa"

#index the assembly
chromap -i -r $assembly -o Atlapetes_index

R1="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Atlapetes_HiC/20260312_LH00541_0124_A22CTYGLT1_SUB17798/fastq/max_hi_c_S1_L001_R1_001.fastq.gz"
R2="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Atlapetes_HiC/20260312_LH00541_0124_A22CTYGLT1_SUB17798/fastq/max_hi_c_S1_L001_R2_001.fastq.gz"

#map the Hi-C reads to the assembly. Output a SAM file
chromap --preset hic -x Atlapetes_index  -r $assembly -1 $R1 -2 $R2 --SAM -o Ablancae.sam

mamba deactivate
mamba activate samtools_env
#Convert SAM file to BAM file
samtools view -bS  Ablancae.sam -o Ablancae.bam

#index assembly
samtools faidx $assembly

#run YaHS:
yahs $assembly Ablancae.bam
