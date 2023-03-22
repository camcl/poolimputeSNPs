#!/bin/bash

# Merge all individually imputed samples into an imputed study population
# Run from appropriate data directory poolimputeSNPs/
# Usage example: apptainer exec container.sif micromamba run -n base bash workflow/scripts/merge_sort_vcf_files.sh results/data/1/prophaser *.full.postgenos.vcf.gz results/data/study.population results/data/1

indir=$1
suffix=$2
samplesId=$3
outdir=$4
outfile=STU.Chr1.SNPs.pruned.sorted.pooled.imputed.vcf.gz

files=$( ls $indir/$suffix )

for file in $files
do
	bcftools index -f $file
done

bcftools merge -Oz -o $outdir/$outfile $files
bcftools view -S $samplesId --force-samples -Oz -o $indir/tmp.sorted.vcf.gz $outdir/$outfile
bcftools view -Oz -o $outdir/$outfile $indir/tmp.sorted.vcf.gz
bcftools index -f $outdir/$outfile
rm $indir/tmp.sorted.vcf.gz 
