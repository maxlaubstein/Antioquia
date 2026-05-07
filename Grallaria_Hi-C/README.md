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
