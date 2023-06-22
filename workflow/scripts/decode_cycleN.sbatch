#!/bin/bash -l
#SBATCH -A snic2022-22-697
# #SBATCH -M snowy
#SBATCH -p core
#SBATCH -C mem128GB
#SBATCH -n 1
#SBATCH -t 00:01:00
#SBATCH -J decodeN
#SBATCH -o slurm-cycle-decodeN-%j.out

# Location: poolimputeSNPs/

# Variables 

ncycle=$1
prev_cycle=$(($ncycle - 1))
max_ncycle=$2
wd=results/data/1

# Run decoding and generate new pooled file

apptainer exec container.sif micromamba run -n base ./opt/simpool/randompoolertune.x cycle$prev_cycle/

