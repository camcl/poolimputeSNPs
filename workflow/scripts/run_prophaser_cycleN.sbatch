#!/bin/bash -l
#SBATCH -A snic2022-22-697
# #SBATCH -M snowy
#SBATCH -p node
#SBATCH -C mem128GB
#SBATCH -n 16
#SBATCH -t 00:01:00
#SBATCH -J prophaserN
#SBATCH -o slurm-cycle-prophaserN-%A_%a.out

# Location: poolimputeSNPs/

# Variables 

ncycle=$1
max_ncycle=$2
wd=results/data/1

echo "Prophaser for sample $SLURM_ARRAY_TASK_ID in cycle $ncycle starting"
echo "Max nb of cycles is $max_ncycle"
echo "*****************************************************************************************************"
echo ""

echo "dummy" >> $wd/dummy-$SLURM_ARRAY_TASK_ID-cycle$ncycle.txt
ls -la $wd
echo ""
df -h $wd
echo ""

echo "dummy" >> $SNIC_TMP/dummy-scratch.txt
ls -la $SNIC_TMP
echo ""
df -h $SNIC_TMP
echo ""

# Impute the sample

apptainer exec container.sif micromamba run -n base bash workflow/scripts/run_workflow_prophaser.sh

# Debugging

status=$?
[ $status -ne 0 ] && echo "apptainer failed with exit code $status"
exit $status
