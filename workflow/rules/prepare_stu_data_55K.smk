import os, sys
from snakemake.utils import min_version
min_version("5.3.0")

"""
Preparation of the data set for the study population.
"""

# TODO: STU_FILE to config

STU_FILE = os.path.join(config["paths"]["data_out"], config["files"][".vcf.gz"]["stu"])  # "results/data/STU.SNPs.pruned.vcf.gz"

# =========================================================================================================
#     Download the data for the inbred lines used as study population from the UCL server
# =========================================================================================================

rule load_stu_data:
	"""
	Download the pruned PLINK files for inbred lines (55K variants) from http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/MAGIC_diverse_FILES/MAGIC_PLINK_PRUNED.tar.gz.
	Create VCF from PLINK files without missing data and chrom 1 -> 21.
	Force major allele to be the reference allele (like  with PLINK 1.x).
	Remove variants with missing genotypes (34495 variants removed).
	"""
	input:
		"opt/genotypooler/genotypooler"
	output:
		plink_gz = "resources/MAGIC_PLINK_PRUNED/ALLchr.SNPs.pruned.vcf.gz",
		plink_vcf = temp("resources/MAGIC_PLINK_PRUNED/ALLchr.SNPs.pruned.vcf")
	params:
		archive_dir = "resources/MAGIC_PLINK_PRUNED",
		plink_suffix = "ALLchr.SNPs.pruned"
	container:
		"plink_1.90b6.21--h516909a_0.sif" # singularity pull "docker://quay.io/biocontainers/plink:1.90b6.21--h516909a_0"
	shell:
		"""
		cd resources
		wget http://mtweb.cs.ucl.ac.uk/mus/www/MAGICdiverse/MAGIC_diverse_FILES/MAGIC_PLINK_PRUNED.tar.gz
	    tar -xvzf MAGIC_PLINK_PRUNED.tar.gz
	    cd MAGIC_PLINK_PRUNED
	    plink --bfile ALLchr.SNPs.pruned --geno 0.0 --recode vcf --out ALLchr.SNPs.pruned
	    cd ../..
	    bcftools view -Oz -o {output.plink_gz} {output.plink_vcf}
	    bcftools index -f {output.plink_gz} 
		"""

# =========================================================================================================
#     Get the coordinates of the target markers in the study population
# =========================================================================================================
	
rule get_coords_vars_stu:
    """
    Extract the physical positions of the 55K variants from PLINK files.
    """
    input:
    	plink_gz = "resources/MAGIC_PLINK_PRUNED/ALLchr.SNPs.pruned.vcf.gz" 
    output:
        os.path.join("results", "data", config["files"][".coords"])
    shell:
        """
        bcftools query -f '%CHROM\t%POS\n' {input.plink_gz} > {output}
        """

# =========================================================================================================
#     Subset the target markers and 496 = 16 * 31 samples for the study population
# =========================================================================================================

rule intersect_coords_samples_stu:
    """
    Intersect coordinates with the entire file and keep only 496 study samples.
    The VCF file used for intersection MUST be indexed (.vcf.gz.csi) and therefore bgzipped (.vcf.gz).
    """
    input:
    	os.path.join("results", "data", config["files"][".coords"])
    output:
        STU_FILE
    params:
        samples_file = os.path.join("results", "data", config["files"][".samples"]),
        archive_dir = "resources/MAGIC_PLINK_PRUNED",
        plink_suffix = "ALLchr.SNPs.pruned"
        # TODO: 496 computed from integer division tot_nb_samples // 16
		# TODO: 16 as block_size param to config
    log:
    	"results/logs/intersect_coords_samples/STU_FILE.log"
    shell:
        """
        bcftools query -l {params.archive_dir}/{params.plink_suffix}.vcf.gz | shuf | head -496 > {params.samples_file}
        bcftools view -R {input} -S {params.samples_file} -Oz -o {output} {params.archive_dir}/{params.plink_suffix}.vcf.gz
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
#     Sort variants by ascending physical position
# =========================================================================================================
        
rule sort_vars_stu:
    """
    Sort variants.
    """
    input:
    	STU_FILE
    output:
        "results/data/STU.SNPs.pruned.sorted.vcf.gz"
    shell:
        """
        bcftools sort -Oz -o {output} {input}
        bcftools index -f {output}
        """

# =========================================================================================================
#     Split data by chromosome
# =========================================================================================================
        
rule split_chrom_stu:
    """
    Split VCF by chromosome.
    """
    input:
    	"results/data/STU.SNPs.pruned.sorted.vcf.gz"
    output:
        "results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.vcf.gz",
        "results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.vcf.gz.csi"
        # NB: this automatically creates the ./{nchrom}/ directory
    params:
        suffix_out = "SNPs.pruned.sorted"
    log:
		"results/logs/split_chrom/{nchrom}.log" 
    shell:
        """
        echo "Writing to results/data/{wildcards.nchrom}/ directory" 
		bcftools view -Oz -o results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz --regions {wildcards.nchrom} {input}
		bcftools index -f results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz
        """
        
# Try to put in log?
'''

		echo "Check dimensions of the data set results/data/{wildcards.nchrom}/Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz and the data fields." > {log}
		echo "Number of samples: " >> {log}
		bcftools query -l results/data/{wildcards.nchrom}/Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz | wc -l >> {log}
		echo "Number of variants: " >> {log}
		bcftools view -H results/data/{wildcards.nchrom}/Chr{wildcards.nchrom}.{params.suffix_out}.vcf.gz | wc -l >> {log}
'''		

