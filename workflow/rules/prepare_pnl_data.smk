import os, sys
from snakemake.utils import min_version
min_version("5.3.0")

"""
Preparation of the data set used as reference panel.
"""

# TODO: PNL_FILE to config

PNL_FILE = "results/data/PNL.SNPs.pruned.vcf.gz"

# =========================================================================================================
#     Download data for the founders from UCL server
# =========================================================================================================

rule load_pnl_data:
	"""
	Download the pruned PLINK files for inbred lines (ca. 1M variants) from http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/MAGIC_diverse_FILES/FOUNDERS.tar.gz.
	Create VCF from PLINK files without missing data and chrom 1 -> 21.
	Force major allele to be the reference allele (like  with PLINK 1.x).
	Remove variants with missing genotypes (0 variant removed).
	"""
	input:
		"opt/genotypooler/genotypooler"
	output:
		plink_gz = "resources/FOUNDERS/Founders.vcf.gz",
		plink_vcf = temp("resources/FOUNDERS/Founders.vcf")
	params:
		archive_dir = "resources/FOUNDERS",
		plink_suffix = "Founders"
	shell:
		"""
		cd resources
		wget http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/MAGIC_diverse_FILES/FOUNDERS.tar.gz
		tar -xvzf FOUNDERS.tar.gz
		cd FOUNDERS
		plink --help
		plink --bfile Founders --geno 0.0 --recode vcf --out Founders
		cd ../..
		bcftools view -Oz -o {output.plink_gz} {output.plink_vcf}
		bcftools index -f {output.plink_gz} 
		"""

# =========================================================================================================
#     Subset from the reference panel the targeted variants in the study population
# =========================================================================================================
		
rule intersect_coords_samples_pnl: 
	"""extract coordinates of the pruned set of markers from the inbred lines file (pruned.nomiss)"""
	input:
		coords_vars = os.path.join("results", "data", config["files"][".coords"]), # input makes sure the file exists
		plink_gz = "resources/FOUNDERS/Founders.vcf.gz" 
	output:
		PNL_FILE
	params: 
		archive_dir = "resources/FOUNDERS",
		plink_suffix = "Founders"
	log:
		"results/logs/intersect_coords_samples/PNL_FILE.log"
	shell:
		"""
		bcftools view -R {input.coords_vars} -Oz -o {output} {input.plink_gz}
		bcftools index -f {output} 
		echo "Check dimensions of the data set {output} and the data fields." > {log}
		echo "Number of samples: " >> {log}
		bcftools query -l {output} | wc -l >> {log}
		echo "Number of variants: " >> {log}
		bcftools view -H {output} | wc -l >> {log}		
		# echo "First variant: " >> {log}
		# bcftools view -H {output} | head -1 | cut -f1-11 >> {log} # this line fails, no idea why
		"""
		
# =========================================================================================================
#     Sort variants
# =========================================================================================================

rule sort_vars_pnl:
	"""
	Sort variants.
	"""
	input:
		PNL_FILE
	output:
		"results/data/PNL.SNPs.pruned.sorted.vcf.gz"
	shell:
		"""
		bcftools sort -Oz -o {output} {input}
		bcftools index -f {output}
		"""

# =========================================================================================================
#     Split data by chromosome
# =========================================================================================================

rule split_chrom_pnl:
	"""
	Split by chromosome.
	"""
	input:
		"results/data/PNL.SNPs.pruned.sorted.vcf.gz"
	output:
		"results/data/{nchrom}/PNL.Chr{nchrom}.SNPs.pruned.sorted.{ext}" 
	params:
		suffix_out = "SNPs.pruned.sorted"
	shell:
		"""
		echo "Writing to results/data/{wildcards.nchrom}/ directory" 
				bcftools view -Oz -o results/data/{wildcards.nchrom}/PNL.Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz --regions {wildcards.nchrom} {input}
				bcftools index -f results/data/{wildcards.nchrom}/PNL.Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz
		"""

