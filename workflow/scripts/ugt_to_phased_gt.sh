#!/bin/bash

# Forces phasing in a VCF file with homozygous unphased GT
# Usage example (from poolimputeSNPs/): 
# $ bash workflow/scripts/ugt_to_phased_gt.sh results/data/tmp.PNL.SNPs.pruned.sorted.vcf.gz results/data/PNL.SNPs.pruned.sorted.vcf.gz

fin=$1
fout=$2
datadir=$( dirname "$fin" )
fbname=$(basename "$fout" .vcf.gz)

sed --help
bcftools view $fin | sed -e 's/1\/1/1\|1/g' -e 's/0\/0/0\|0/g' -e 's/\.\/\./\.\|\./g' > $datadir/$fbname.vcf
bcftools view -Oz -o $fout $datadir/$fbname.vcf
bcftools index -f $fout
rm $datadir/$fbname.vcf

