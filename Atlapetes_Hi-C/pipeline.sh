mamba activate chromap_env

assembly="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa"

#Teresa used fastp to remove duplicates prior to mapping

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

#add read groups for picard
mamba activate picard_env
picard AddOrReplaceReadGroups I=Ablancae.bam O=Ablancae_rg.bam RGID=1 RGLB=HiC RGPL=ILLUMINA  RGPU=unit1 RGSM=Ablancae

#sort BAM file
mamba activate samtools_env
samtools sort -@ 16 -o Ablancae_rg_sorted.bam Ablancae_rg.bam

#sort again by first column
samtools sort -n -@ 16 -o Ablancae_YaHS.bam Ablancae_marked_dups.bam

#index assembly
samtools faidx $assembly
samtools faidx /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/blancae.p_ctg.fa

#run YaHS:
yahs $assembly Ablancae_YaHS.bam
