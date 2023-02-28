#!/bin/bash


# all paths are relative phaser/   (assumed this script is executed from workflow/scripts)
exe_path=../../opt/prophaser
cd $exe_path

ne=16
error=1e-12
mapfile=../../results/data/1/1_iwgsc_refseqv1.0_mapping_data
indir=../../results/data/1/
sample=$SLURM_ARRAY_TASK_ID
samples_file=../../results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz
ref_file=../../results/data/1/PNL.Chr1.SNPs.pruned.sorted.vcf.gz
results_directory=../../results/data/1/prophaser/

# parameters for the linear_state  branch
algo=integrated
algo=separate
algo=fixed
niter=3

# write 1 file for every sample
sampleId=$( bcftools query -l $indir$samples_file | head -$sample | tail -1 )
sample_file=s$sampleId.$samples_file 
bcftools view -Oz -s $sampleId -o $indir$sample_file $indir$samples_file


# prepare directory and template
if [ ! -d "$results_directory" ]; then
mkdir -p $results_directory
fi

./create_template_vcf.sh $indir $sample_file
./create_template_vcf_gtgp.sh $indir $sample_file


export OMP_NUM_THREADS=16 # Snowy

##### assuming the master branch has been compiled to this executable
./phase --results_directory $results_directory  --directory $indir --sample_file $sample_file --reference_file $ref_file --ne $ne --map $mapfile --error $error 


# clean file
rm $indir$sample_file
