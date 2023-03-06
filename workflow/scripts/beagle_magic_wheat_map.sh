#!/bin/bash

PNL=results/data/1/PNL.Chr1.SNPs.pruned.sorted
STU=results/data/1/STU.Chr1.SNPs.pruned.sorted.pooled
map=results/data/1/1_interpolated_wheat_map  
# 1_interpolated_wheat_map_plink.map
outdir=results/data/1/beagle

beaglejar=workflow/scripts/beagle.11Mar19.69c.jar
cfgtjar=workflow/scripts/conform-gt.jar

modelscale=0.8
ne=1000000
niterations=5

chrom=$( bcftools query -f '%CHROM\n' $PNL.vcf.gz | head -1 )
startpos=$( bcftools query -f '%POS\n' $PNL.vcf.gz | head -1 )
endpos=$( bcftools query -f '%POS\n' $PNL.vcf.gz | tail -1 )

if [ ! -d "$outdir" ]; then
mkdir -p $outdir
fi

echo 'Contigs in the reference file'
echo '.................................................................................'
echo 'Chromosome ' $chrom '   Startpos =' $startpos '   Endpos =' $endpos 
echo ''

echo ''
echo 'Check FORMAT field in files for imputation'
echo '.................................................................................'
reftruefmt=$( bcftools query -f '%LINE\n' $PNL.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in reference panel: ' $reftruefmt
imppoolfmt=$( bcftools query -f '%LINE\n' $STU.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in target: ' $imppoolfmt
echo ''

echo ''
echo 'Check number of samples and number of markers in files for imputation'
echo '.................................................................................'
echo 'reference:'
bcftools query -l $PNL.vcf.gz | wc -l
#bcftools view -H  $PNL.vcf.gz | wc -l
echo ''
echo 'target:'
bcftools query -l $STU.vcf.gz | wc -l
#bcftools view -H $STU.vcf.gz | wc -l
echo ''

echo ''
echo 'Phase reference and target with BEAGLE'
echo '.................................................................................'
echo 'Beagle .jar file used at:' $beaglejar
echo ''
#java -Xss5m -jar $beaglejar impute=false gtgl=$PNL.vcf.gz out=$PNL.phased
cp $PNL.vcf.gz $PNL.phased.vcf.gz
bcftools index -f $PNL.phased.vcf.gz
refphasfmt=$( bcftools query -f '%LINE\n' $PNL.phased.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the phased ref file:' $refphasfmt 

### two runs for phasing from GL: 1st run from GL to GT:DP:GP (unphased GT)
java -Xss5m -jar $beaglejar impute=false gtgl=$STU.vcf.gz map=$map modelscale=$modelscale niterations=$niterations ne=$ne out=$STU.unphased
### two runs for phasing from GL: 2nd run from GT:DP:GP to phased GT
### with gt argument, all genotypes in the output file will be phased and non-missing  (Beagle4.1 documentation)
echo ''
java -Xss5m -jar $beaglejar impute=false gt=$STU.unphased.vcf.gz map=$map modelscale=$modelscale niterations=$niterations ne=$ne out=$STU.phased
bcftools index -f $STU.phased.vcf.gz
impphasfmt=$( bcftools query -f '%LINE\n' $STU.phased.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the phased target file:' $impphasfmt
echo ''

echo ''
echo 'Deduplicate possibly duplicated markers'
echo '.................................................................................'
bcftools norm --rm-dup all -Oz -o $STU.phased.dedup.vcf.gz $STU.phased.vcf.gz
bcftools index -f $STU.phased.dedup.vcf.gz
bcftools norm --rm-dup all -Oz -o $PNL.phased.dedup.vcf.gz $PNL.phased.vcf.gz
bcftools index -f $PNL.phased.dedup.vcf.gz
echo ''

echo ''
echo 'Unify reference and target markers with CONFORM-GT'
echo '.................................................................................'
echo 'conform-gt .jar file used at:' $cfgtjar
### necessary to get proper imputation results and GT:DS:GP fields with gprobs=true
echo ''
java -jar $cfgtjar ref=$PNL.phased.dedup.vcf.gz gt=$STU.phased.dedup.vcf.gz chrom=$chrom:$startpos-$endpos out=$STU.cfgt
bcftools index -f $STU.cfgt.vcf.gz
echo ''

echo ''
echo 'Impute target from reference with BEAGLE'
echo '.................................................................................'
echo 'Beagle .jar file used at:' $beaglejar
### impute=true must be used with gt= and ref= (Beagle4.1 documentation)
echo ''
### map=$chrom_interpolated_wheat_map_plink.map does not work (why???), looks like it has to be hard cpded
java -Xss5m -jar $beaglejar gt=$STU.cfgt.vcf.gz ref=$REF.phased.vcf.gz map=$map modelscale=$modelscale niterations=$niterations ne=$ne impute=true gprobs=true out=$STU.imputed
bcftools index -f $STU.imputed.vcf.gz
impimpfmt=$( bcftools query -f '%LINE\n' $STU.imputed.vcf.gz | head -1 | cut -f9 )
echo 'FORMAT in the imputed target file:' $impimpfmt
cp $STU.imputed.vcf.gz $outdir/
cp $STU.imputed.vcf.gz.csi $outdir/
echo''

echo ''
echo 'Cleaning directory from log files'
echo '.................................................................................'
rm ./*.log
rm ./*cfgt.vcf.gz* # files from conform-gt cannot be overwritten
rm ./*phased*.vcf.gz*
rm ./*imputed.vcf.gz*
echo 'Done.'
echo ''
