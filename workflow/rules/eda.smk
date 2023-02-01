import os, sys

from snakemake.utils import min_version
min_version("5.3.0")

# =========================================================================================================
#     Remove markers with swapped alleles
# =========================================================================================================

rule remove_swaps:
	"""
	Remove marker positions where the alternate allele is swapped between
	the reference and the study populations.
	
	Mess in the packages and modules here!
	Necessary to install with `pip`all the requrements for genotypooler:
	```
	(smkenv) camille@camille-Precision-3551:~/MagicWheat/runs/devtests$ python3 -V
Python 3.6.15
(smkenv) camille@camille-Precision-3551:~/MagicWheat/runs/devtests$ pip3 install -r /home/camille/MagicWheat/src/genotypooler/requirements.txt

	```
	"""
	input:
		"results/data/{nchrom}/PNL.Chr{nchrom}.SNPs.pruned.sorted.vcf.gz",
		"results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.vcf.gz",
		"opt/genotypooler/genotypooler/runtools/rm_swapped_ref_alt.py"
		
	output:
		"results/data/{nchrom}/unfiltered-opposite-markers/swapped_chrom_pos.coords"
		
	params:
		path_to_script = "opt/genotypooler/genotypooler/runtools",
		cwd = os.getcwd()
		
	log:
		os.path.join(os.getcwd(),  "results/logs/remove_swaps/{nchrom}.log")  # python outputs
		
	shell:
		"""
		cd {params.path_to_script}
		python3 -u rm_swapped_ref_alt.py {params.cwd}/results/data/{wildcards.nchrom} {params.cwd}/results/data/{wildcards.nchrom} PNL.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz > {log}
		cd {params.cwd}
		"""

# =========================================================================================================
#     Plot histograms of marker composition
# =========================================================================================================

rule plot_histograms:
	"""
	Plot histograms showing the composition in markers of the data sets for the study population and the reference panel.
	"""
	input:
		"results/data/{nchrom}/unfiltered-opposite-markers/swapped_chrom_pos.coords"
		
	output:
		"results/plots/{nchrom}/PNL/genotypes_hexa_scaled_proportions.pdf",
		"results/plots/{nchrom}/STU/genotypes_hexa_scaled_proportions.pdf"
		
	params:
		path_to_script = 'opt/genotypooler/genotypooler/graphtools',
		cwd = os.getcwd(),
		dir_out_pnl = directory("results/plots/{nchrom}/PNL"),
		dir_out_stu = directory("results/plots/{nchrom}/STU")
		
	log:
		os.path.join(os.getcwd(),  "results/logs/plot_histograms/{nchrom}.log")  # python outputs
		
	shell:
		"""
		cd {params.path_to_script}
		python3 -u genotypes_decoding.py {params.cwd}/results/data/{wildcards.nchrom}/PNL.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/results/data/{wildcards.nchrom}/PNL.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/{params.dir_out_pnl}  > {log}
		python3 -u genotypes_decoding.py {params.cwd}/results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/{params.dir_out_stu} > {log}
		cd {params.cwd}
		"""
