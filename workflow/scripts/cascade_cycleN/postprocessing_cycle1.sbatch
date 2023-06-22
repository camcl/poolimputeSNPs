#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -C mem128GB
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J postcycle1
#SBATCH -o slurm-cycle-post1-%j.out

# Location: poolimputeSNPs/

# Variables

ncycle=1
wd=results/data/1
LIBSTATGEN=$( apptainer exec container.sif micromamba run -n base ls /opt/conda/pkgs | grep libstatgen )
echo "LIBSTATGEN=$LIBSTATGEN"
echo ""

# Merge imputed data 

apptainer exec container.sif micromamba run -n base bash workflow/scripts/merge_sort_vcf_files.sh results/data/1/prophaser *.full.postgenos.vcf.gz results/data/study.population results/data/1

# Clean directory

mkdir $wd/cycle$ncycle
mv $wd/prophaser $wd/cycle$ncycle/ && mv $wd/*pooled* $wd/cycle$ncycle/

# Create GL/GP and log-GL files from true GT (necessary for randompoolertune and for later plotting)

bash workflow/scripts/ugt_to_log_gl.sh ./results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz ./results/data/1/STU.Chr1.SNPs.pruned.sorted.log-gl.vcf.gz
bash workflow/scripts/ugt_to_gl.sh ./results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz ./results/data/1/STU.Chr1.SNPs.pruned.sorted.gl.vcf.gz

# Compile randompoolertune

apptainer exec container.sif micromamba run -n base g++ opt/simpool/randompoolertune.cpp -O3 -I /opt/conda/pkgs/$LIBSTATGEN/include -lStatGen -L /opt/conda/pkgs/$LIBSTATGEN/ -lz -o opt/simpool/randompoolertune.x

# Ready for starting next cycle

echo ""
echo "Postprocessing of cycle 1 done"
echo "*****************************************************************************************************"
echo ""
sbatch --dependency=afterok:$SLURM_JOB_ID workflow/scripts/decode_cycleN.sbatch 2 15
