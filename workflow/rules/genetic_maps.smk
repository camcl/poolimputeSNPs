import os
import pandas as pd
from snakemake.utils import min_version
min_version("5.3.0")


# =========================================================================================================
#     After pooling: create the genetic maps required for the post-pooling imputation step
# =========================================================================================================

rule load_mapping_data:
	"""
	Download and unzip the data from URGI INRA server: 
	https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_recombination_rate_analysis.zip.
	The map was created based on the RefSeq v1.0 genome assembly by IWGSC,
	which is the assembly that was used for aligning the reads in the NIAB Diverse MAGIC wheat population. 
	"""
	input:		
		eda_pnl_figs = expand("results/plots/{nchrom}/PNL/genotypes_hexa_scaled_proportions.pdf", nchrom=list(config["chromosomes"]["prefix"].values())), 
		eda_stu_figs = expand("results/plots/{nchrom}/STU/genotypes_hexa_scaled_proportions.pdf", nchrom=list(config["chromosomes"]["prefix"].values())) 
	output:
		directory("resources/iwgsc_refseqv1.0_recombination_rate_analysis")
	log:
		os.path.join(os.getcwd(), "results/logs/load_genetic_map")
	shell:
		"""
		cd resources
		wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_recombination_rate_analysis.zip > {log}
		wget https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_recombination_rate_analysis.zip.md5.txt >> {log}
		unzip iwgsc_refseqv1.0_recombination_rate_analysis.zip >> {log}
		"""


rule split_map_by_chromosome: 
	"""
	Create one map per chromosome (required for imputation).
	Note that some chromosomes might not have any marker mapped, in such case an empty file is created.
	"""
	input:
		recomb_map = "resources/iwgsc_refseqv1.0_recombination_rate_analysis"
	output:
		chrom_recomb_map = "results/data/{nchrom}/{nchrom}_iwgsc_refseqv1.0_mapping_data"
	params:
		map_file = "iwgsc_refseqv1.0_mapping_data.txt",
		bwd_chrom_tab = config["chromosomes"]["bwd"]
	log:
		write_to = os.path.join(os.getcwd(), "results/logs/split_map_by_chromosome/{nchrom}.log")
	run:
		print(f'Working on chromosome {wildcards.nchrom} = chr{params.bwd_chrom_tab[int(wildcards.nchrom)]}'.ljust(80, '.'))
		all_chr_data = pd.read_csv(os.path.join(input.recomb_map, params.map_file), sep='\t')
		map_data = all_chr_data[all_chr_data['chromosome'] == f'chr{params.bwd_chrom_tab[int(wildcards.nchrom)]}']
		map_data.to_csv(output.chrom_recomb_map,
				    sep='\t',
				    columns=['psId', 'physicalPosition', 'geneticPosition'],
				    header=False,
				    index=False)
		with open(log.write_to, 'w') as f_log:
			print(map_data, file=f_log)
			
			
rule get_marker_id_pos:
	"""
	Write for every chromosome to a text file the ID and the POS of each marker (1 marker per line) in the reference panel.
	"""
	input:
		chrom_recomb_map = "results/data/{nchrom}/{nchrom}_iwgsc_refseqv1.0_mapping_data",
		pnl_vcf = "results/data/{nchrom}/PNL.Chr{nchrom}.SNPs.pruned.sorted.vcf.gz"
	output:
		file_id_pos = "results/data/{nchrom}/{nchrom}_markers_is_pos.txt"
	shell:
		"""
		bcftools query -f "%ID\t%POS\n" {input.pnl_vcf} > {output.file_id_pos}
		"""
		

rule interpolate_chrom_map: 
	"""
	Perform for each chromosome an interpolation of the genetic distances for the unmapped markers from the reference panel.
	The interpolated genetic position must have ascending values for compatibility with Beagle4.1.
	In the case no marker was mapped on the chromosome, the genetic distances between markers are set to be 1 cM/Mb as suggested for Beagle 4.1.
	"""
	input:
		file_id_pos = "results/data/{nchrom}/{nchrom}_markers_is_pos.txt",
		chrom_recomb_map = "results/data/{nchrom}/{nchrom}_iwgsc_refseqv1.0_mapping_data"
	output:
		chrom_interpol_map = "results/data/{nchrom}/{nchrom}_interpolated_wheat_map",
		tmp_map = temp("results/data/{nchrom}/{nchrom}_interpolated_wheat_map_tmp")
	params:
		map_file = "iwgsc_refseqv1.0_mapping_data.txt",
		bwd_chrom_tab = config["chromosomes"]["bwd"]
	log:
		write_to = os.path.join(os.getcwd(), "results/logs/interpolate_chrom_map/{nchrom}.log")
	script:
		"../scripts/interpolate_maps.py"

