#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -C mem128GB
#SBATCH -n 1
#SBATCH -t 00:01:00
#SBATCH -J imputeN
#SBATCH -o slurm-cycle-imputeN-%j.out

# Location: poolimputeSNPs/

# Variables 

ncycle=$1
max_ncycle=$2
wd=results/data/1

echo "Imputation in cycle $ncycle starting"
echo "Max nb of cycles is $max_ncycle"
echo "*****************************************************************************************************"
echo ""

# Deploy array job for imputation

echo ""
echo "*****************************************************************************************************"
echo ""
sbatch --array=1-496%1 workflow/scripts/run_prophaser_cycleN.sbatch $ncycle $max_ncycle
