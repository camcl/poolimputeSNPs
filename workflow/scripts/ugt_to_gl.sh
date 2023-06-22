#!/bin/bash

# Converts a VCF file written with unphased GT to a VCF file with log GL
# Usage example (from poolimputeSNPs/): 
# $ bash ../ugt_to_gl.sh /crex/proj/snic2019-8-216/private/iter-pool-imp/runs/poolimputeSNPs/results/data/1/STU.Chr1.SNPs.pruned.sorted.vcf.gz /crex/proj/snic2019-8-216/private/iter-pool-imp/runs/poolimputeSNPs/results/data/1/STU.Chr1.SNPs.pruned.sorted.gl.vcf.gz

# Updates the metadata in the header
# $ ~/PoolImpHuman/data/20200615$ bcftools view -h IMP.chr20.snps.gt.chunk10000.vcf.gz | sed 's/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">/##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">/' | grep '##FORMAT'
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">

# Replaces GT FORMAT field with GL for each variant
# $ bcftools view -H IMP.chr20.snps.gt.chunk10000.vcf.gz | sed 's/GT/GL/g'

# Converts the GT values to LOG GL ones, including missing
# $ bcftools view -H IMP.chr20.snps.gt.chunk10000.vcf.gz | sed -e ...

fin=$1
fout=$2
datadir=$( dirname "$fin" )
fbname=$(basename "$fout" .vcf.gz)

bcftools view $fin | sed 's/##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">/##FORMAT=<ID=GL,Number=G,Type=Float,Description="Estimated Genotype Probability">/' | sed 's/GT/GL/g' | sed -e 's/1\/1/0,0,1/g' -e 's/0\/0/1,0,0/g' -e 's/1\/0/0,1,0/g' -e 's/0\/1/0,1,0/g' -e 's/.\/./0.333333333,0.333333333,0.333333333/g' > $datadir/$fbname.vcf
bcftools view -Oz -o $fout $datadir/$fbname.vcf
bcftools index -f $fout
rm $datadir/$fbname.vcf
