import os, sys

from snakemake.utils import min_version
min_version("5.3.0")

# =========================================================================================================
#     Plot histograms of marker composition
# =========================================================================================================

rule plot_histograms:
	"""
	Plot histograms showing the composition in markers of the data sets for the study population and the reference panel.
	"""
	input:
		expand("results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.pooled.vcf.gz", nchrom=list(config["chromosomes"]["prefix"].values())),
		expand("results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.pooled.vcf.gz.csi", nchrom=list(config["chromosomes"]["prefix"].values()))
		
	output:
		"results/plots/{nchrom}/PNL/genotypes_hexa_scaled_proportions.pdf",
		"results/plots/{nchrom}/STU/genotypes_hexa_scaled_proportions.pdf",
		"results/plots/{nchrom}/pooledSTU/genotypes_hexa_scaled_proportions.pdf"
		
	params:
		path_to_script = 'opt/genotypooler/genotypooler/graphtools',
		cwd = os.getcwd(),
		dir_out_pnl = directory("results/plots/{nchrom}/PNL"),
		dir_out_stu = directory("results/plots/{nchrom}/STU"),
		dir_out_pool = directory("results/plots/{nchrom}/pooledSTU")
		
	log:
		os.path.join(os.getcwd(),  "results/logs/plot_histograms/{nchrom}.log")  # python outputs
		
	shell:
		"""
		cd {params.path_to_script}
		python3 -u genotypes_distrib_hist.py {params.cwd}/results/data/{wildcards.nchrom}/PNL.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/results/data/{wildcards.nchrom}/PNL.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/{params.dir_out_pnl}  > {log}
		python3 -u genotypes_distrib_hist.py {params.cwd}/results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/{params.dir_out_stu} >> {log}
		python3 -u genotypes_distrib_hist.py {params.cwd}/results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {params.cwd}/results/data/{wildcards.nchrom}/STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.pooled.vcf.gz {params.cwd}/{params.dir_out_pool} >> {log}
		cd {params.cwd}
		"""
