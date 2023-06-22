#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -C mem128GB
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -J Ncycles
#SBATCH -o slurm-cycleN-run-%j.out

# Location: poolimputeSNPs/
# If this master job is submitted on Snowy, all other depending job downstream are run on Snowy as well.

# Variables

ncycle=41
max_ncycle=55
wd=results/data/1
LIBSTATGEN=$( apptainer exec container.sif micromamba run -n base ls /opt/conda/pkgs | grep libstatgen )

# Create GL/GP and log-GL files from true GT (necessary for randompoolertune and for later plotting)
echo ""
echo "Converting true data to GL/GP and log-GL"
echo "*****************************************************************************************************"
echo ""

bash workflow/scripts/ugt_to_log_gl.sh ./results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz ./results/data/1/STU.Chr1.SNPs.pruned.sorted.log-gl.vcf.gz
bash workflow/scripts/ugt_to_gl.sh ./results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz ./results/data/1/STU.Chr1.SNPs.pruned.sorted.gl.vcf.gz


# Compile randompoolertune

echo "LIBSTATGEN=$LIBSTATGEN"
echo ""
echo "Compiling randompoolertune"
echo "*****************************************************************************************************"
echo ""

apptainer exec container.sif micromamba run -n base g++ opt/simpool/randompoolertune.cpp -O3 -I /opt/conda/pkgs/$LIBSTATGEN/include -lStatGen -L /opt/conda/pkgs/$LIBSTATGEN/ -lz -o opt/simpool/randompoolertune.x

jid_prev=$(sbatch --parsable workflow/scripts/postprocessing_cycleN.sbatch 1 $max_ncycle)

# Ready for starting next cycles

for N in {41..55}
do
	echo ""
	echo "Decoding in cycle $N/$max_ncycle scheduled"
	echo "*****************************************************************************************************"
	echo ""
	jid_decodeN=$(sbatch --parsable --dependency=afterok:${jid_prev} workflow/scripts/decode_cycleN.sbatch $N $max_ncycle)
	echo ""
	
	echo "Imputation in cycle $N/$max_ncycle scheduled"
	echo "*****************************************************************************************************"
	echo ""
	jid_imputeN=$(sbatch --parsable --dependency=afterok:${jid_prev}:${jid_decodeN} --array=1-496 workflow/scripts/run_prophaser_cycleN.sbatch $N $max_ncycle)
	echo ""

	echo "Postprocessing in cycle $N/$max_ncycle scheduled"
	echo "*****************************************************************************************************"
	echo ""
	jid_postN=$(sbatch --parsable --dependency=afterok:${jid_prev}:${jid_decodeN}:${jid_imputeN} workflow/scripts/postprocessing_cycleN.sbatch $N $max_ncycle)
	echo ""

	jid_prev=$jid_postN
done




