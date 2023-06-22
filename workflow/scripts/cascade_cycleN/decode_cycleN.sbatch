#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -C mem128GB
#SBATCH -n 1
#SBATCH -t 00:02:00
#SBATCH -J decodeN
#SBATCH -o slurm-cycle-decodeN-%j.out

# Location: poolimputeSNPs/

# Variables 

ncycle=$1
prev_cycle=$(($ncycle - 1))
max_ncycle=$2
wd=results/data/1

echo "Cycle $ncycle starting"
echo "Max nb of cycles is $max_ncycle"
echo "*****************************************************************************************************"
echo ""

# Run decoding and generate new pooled file

if [ $ncycle -le $max_ncycle ]; then
	apptainer exec container.sif micromamba run -n base ./opt/simpool/randompoolertune.x cycle$prev_cycle/
	echo "*****************************************************************************************************"
	echo ""
	sbatch --dependency=afterok:$SLURM_JOB_ID workflow/scripts/impute_cycleN.sbatch $ncycle $max_ncycle
else 
	echo ""
	echo "Max number of cycles is reached, exiting now"
	echo "*****************************************************************************************************"
	echo ""
	exit
fi
