#!/bin/bash

# Merge all individually imputed samples into an imputed study population
# Run from appropriate data directory poolimputeSNPs/results/data/1
# Usage example: bash ../../../workflow/scripts/merge_sort_vcf_files.sh ./prophaser *.full.postgenos.vcf.gz ../study.population

indir=$1
suffix=$2
samplesId=$3
outdir=$( pwd )
outfile=STU.Chr1.SNPs.pruned.sorted.pooled.imputed.vcf.gz

module load bioinfo-tools && module load bcftools/1.9

files=$( ls $indir/$suffix )

for file in $files
do
	bcftools index -f $file
done

bcftools merge -Oz -o $outdir/$outfile $files
bcftools view -S $samplesId -Oz -o tmp.sorted.vcf.gz $outfile
bcftools view -Oz -o $infile tmp.sorted.vcf.gz
bcftools index $outdir/$outfile
rm tmp.sorted.vcf.gz 
