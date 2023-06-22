#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -C mem128GB
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -J postN
#SBATCH -o slurm-cycle-postN-%j.out

# Location: poolimputeSNPs/

# Variables 

ncycle=$1
max_ncycle=$2
wd=results/data/1

echo "Postprocessing in cycle $ncycle starting"
echo "Max nb of cycles is $max_ncycle"
echo "*****************************************************************************************************"
echo ""

# Merge imputed data 

apptainer exec container.sif micromamba run -n base bash workflow/scripts/merge_sort_vcf_files.sh results/data/1/prophaser *.full.postgenos.vcf.gz results/data/study.population results/data/1

# Clean directory

mkdir $wd/cycle$ncycle
mv $wd/prophaser $wd/cycle$ncycle/
mv $wd/*pooled* $wd/cycle$ncycle/

# Ready for starting next cycle

echo ""
echo "Cycle $ncycle done"
echo "*****************************************************************************************************"
echo ""

sbatch --dependency=afterok:$SLURM_JOB_ID workflow/scripts/decode_cycleN.sbatch $(($ncycle + 1)) $max_ncycle


