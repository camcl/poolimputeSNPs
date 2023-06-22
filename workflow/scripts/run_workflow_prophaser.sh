#!/bin/bash


# all paths are relative phaser/   (assumed this script is executed from poolimputeSNPs/)
proj_root=$(pwd)
exe_path=opt/prophaser
cd $exe_path

echo "Parsing parameters for prophaser"
echo "**********************************************************************************"
ne=16
error=1e-12
mapfile=$proj_root/results/data/1/1_interpolated_wheat_map
indir=$proj_root/results/data/1/
sample=$SLURM_ARRAY_TASK_ID
samples_file=STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz
ref_file=$proj_root/results/data/1/PNL.Chr1.SNPs.pruned.sorted.vcf.gz
results_directory=$proj_root/results/data/1/prophaser/

# write 1 file for every sample
sampleId=$( bcftools query -l $indir$samples_file | head -$sample | tail -1 )
sample_file=s$sampleId.$samples_file 
bcftools view -Oz -s $sampleId -o $indir$sample_file $indir$samples_file

echo ""
echo "Sample processed: $sampleId (nr. $sample)"
echo ""

# prepare directory and template
if [ ! -d "$results_directory" ]; then
mkdir -p $results_directory
fi

./create_template_vcf.sh $indir $sample_file
./create_template_vcf_gtgp.sh $indir $sample_file


export OMP_NUM_THREADS=16 # Snowy

##### assuming the multilevel branch has been compiled to this executable
./phase --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file --ne $ne --map $mapfile --error $error 


# clean file
rm $indir$sample_file
rm $indir"s"$sampleId.*.template*.vcf
