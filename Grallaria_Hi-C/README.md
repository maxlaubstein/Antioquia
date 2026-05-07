# Grallaria Hi-C

## Index the assembly with BWA:

```cd /n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Arima_HiC_Pipeline_Grallaria```

```index.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=index
#SBATCH --partition=shared,edwards
#SBATCH --time=10:00:00
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

## Also index with samtools:

```faidx.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=faidx
#SBATCH --partition=shared,edwards
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/faidx.%j.log
#SBATCH --error=logs/faidx.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate samtools_env

mkdir -p assembly
mkdir -p logs
assembly="assembly/grallaria.p_ctg.fa"
samtools faidx $assembly
~~~


## Trim the Hi-C reads:

```trim.sbatch```:
~~~
 #!/bin/bash
#SBATCH --job-name=trim
#SBATCH --partition=shared,edwards
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/trim.%j.log
#SBATCH --error=logs/trim.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate cutadapt

mkdir -p trimmed_reads
#trim adapters off R1 and R2:
R1="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/HiC_Reads/100bp/20260312_LH00541_0124_A22CTYGLT1_SUB17798/fastq/max_hi_c_S1_L001_R1_001.fastq.gz"
R2="/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/HiC_Reads/100bp/20260312_LH00541_0124_A22CTYGLT1_SUB17798/fastq/max_hi_c_S1_L001_R2_001.fastq.gz"

cutadapt -j 16 -u 5 -U 5 -o trimmed_reads/R1.trimmed.fastq.gz -p trimmed_reads/R2.trimmed.fastq.gz $R1 $R2
~~~

## Map the R1 reads:

```map_R1.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=map_R1
#SBATCH --partition=shared,edwards
#SBATCH --time=10:00:00
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

## Map the R2 reads:

```map_R2.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=map_R2
#SBATCH --partition=shared,edwards
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/map_R2.%j.log
#SBATCH --error=logs/map_R2.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate samtools_env

assembly="assembly/grallaria.p_ctg.fa"
#trimmed read paths:
R2_trimmed="trimmed_reads/R2.trimmed.fastq.gz"

mkdir -p aligned_bam
~/bwa/bwa mem -t 16 $assembly $R2_trimmed | samtools view -@ 16  -Sb - > aligned_bam/grallaria_R2.bam

#get stats
samtools flagstat aligned_bam/grallaria_R2.bam > grallaria_R2.bam.stats
~~~

Now it's time for more work. In this directory ```/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/Antioquia/Arima_HiC_Pipeline_Grallaria``` I have another directory ```Arima_Scripts``` containing several scripts provided by the Arima Mapping Pipeline (https://github.com/ArimaGenomics/mapping_pipeline/blob/master/filter_five_end.pl). I literally just copy-pasted the scripts in.

## Filter 5' ends of chimeric reads (R1):

```filter_five_R1.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=filter_5_R1
#SBATCH --partition=shared,edwards
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/filter_5_R1.%j.log
#SBATCH --error=logs/filter_5_R1.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate samtools_env

mkdir -p filtered_bam

samtools view -h  aligned_bam/grallaria_R1.bam  | perl Arima_Scripts/filter_five_end.pl | samtools view -Sb - > filtered_bam/grallaria_R1.bam
~~~

## Filter 5' ends of chimeric reads (R2):

```filter_five_R2.sbatch```:
~~~
#!/bin/bash
#SBATCH --job-name=filter_5_R2
#SBATCH --partition=shared,edwards
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/filter_5_R2.%j.log
#SBATCH --error=logs/filter_5_R2.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

source ~/.bashrc
mamba activate samtools_env

mkdir -p filtered_bam

samtools view -h  aligned_bam/grallaria_R2.bam  | perl Arima_Scripts/filter_five_end.pl | samtools view -Sb - > filtered_bam/grallaria_R2.bam
~~~

## Pair single-end read reads:

```pair_single_ends.sbatch```:
~~~
~~~
