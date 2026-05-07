# Grallaria Hi-C

## Index the assembly:

```cd /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Arima_HiC_Pipeline_Grallaria```

```index.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --partition=shared,edwards
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/index.%j.log
#SBATCH --error=logs/index.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

mkdir -p assembly
mkdir -p logs
cp /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/genomes/assembly/grallaria.p_ctg.fa assembly/
assembly="assembly/grallaria.p_ctg.fa"

#index assembly:
~/bwa/bwa index -a bwtsw $assembly
~~~

## Trim the Hi-C reads:

```trim.sbatch```:
~~~
 #!/bin/bash
#SBATCH --job-name=index
#SBATCH --partition=shared,edwards
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/trim.%j.log
#SBATCH --error=logs/trim.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate cutadapt

mkdir -p trimmed_reads
#trim adapters off R1 and R2:
R1="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Atlapetes_HiC/20260312_LH00541_0124_A22CTYGLT1_SUB17798/fastq/max_hi_c_S1_L001_R1_001.fastq.gz"
R2="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Atlapetes_HiC/20260312_LH00541_0124_A22CTYGLT1_SUB17798/fastq/max_hi_c_S1_L001_R2_001.fastq.gz"

cutadapt -j 16 -u 5 -U 5 -o trimmed_reads/R1.trimmed.fastq.gz -p trimmed_reads/R2.trimmed.fastq.gz $R1 $R2

~~~

## Map the R1 reads:

```map_R1.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=map_R1
#SBATCH --partition=shared,edwards
#SBATCH --time=60:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/map_R1.%j.log
#SBATCH --error=logs/map_R1.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate samtools_env

assembly="assembly/grallaria.p_ctg.fa"
#trimmed read paths:
R1_trimmed="trimmed_reads/R1.trimmed.fastq.gz"

mkdir -p aligned_bam
~/bwa/bwa mem -t 16 $assembly $R1_trimmed | samtools view -@ 16  -Sb - > aligned_bam/grallaria_R1.bam

#get stats
samtools flagstat aligned_bam/grallaria_R1.bam > grallaria_R1.bam.stats




~~~
