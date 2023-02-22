import os, sys

rootdir = 'opt/genotypooler'
sys.path.insert(0, rootdir)

from snakemake.utils import min_version
min_version("5.3.0")


# =========================================================================================================
#     Chunk chromosomes into subfiles of 1000 markers
# =========================================================================================================
	
rule chunk_chromosomes: 
	"""Pack the study population data for each chromosome into chunks of consecutive markers."""
	# TODO: add verif overall nb markers per chunk sum up to tot nb markers on chromosome
	# TODO: chunk_size to config
	input: # output from EDA
		eda_pnl_figs = "results/plots/{nchrom}/PNL/genotypes_hexa_scaled_proportions.pdf", 
		eda_stu_figs = "results/plots/{nchrom}/STU/genotypes_hexa_scaled_proportions.pdf" 
	output:
		directory("results/data/{nchrom}/tmp")
	params:
		chunk_size = 1000
	log:
		os.path.join(os.getcwd(), "results/logs/chunk_chromosomes/{nchrom}.log")
	shell:
		'''
		cd results/data/{wildcards.nchrom}
		echo "Chunk data for chromosome {wildcards.nchrom}..." > {log}
		echo "" >> {log}
		bash ../../../opt/genotypooler/bin/bcfchunkpara.sh STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz ./tmp {params.chunk_size} >> {log}
		echo "" >> {log}
		echo "Verify chunked files:" >> {log}
		for vcfgz in ./tmp/*.vcf.gz
		do
			echo "Number of markers in $vcfgz" >> {log}
			bcftools view -H $vcfgz | wc -l >> {log}
			echo "Number of samples in $vcfgz" >> {log}
			bcftools query -l $vcfgz | wc -l >> {log}
			echo "" >> {log}
		done
		echo "Leaving directory {wildcards.nchrom}." >> {log}
		cd ../..
		'''

# =========================================================================================================
#     Simulate pooling for each chunked file and merge pooled chunks by chromosome
# =========================================================================================================
	
rule pool_chromosomes:
	"""Apply pooling simulation to the per-chromosome marker data in the study population."""
	# TODO: nb_cores as hyperparameter (nb_threads?)
	input: 
		"results/data/{nchrom}/tmp",
		"resources/adaptive_gls.csv"
	output: 
		"results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.pooled.vcf.gz",
		"results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.pooled.vcf.gz.csi"
	params:
		genotypooler = 'opt/genotypooler/genotypooler',
		cwd = os.getcwd(),
		nb_cores = 4
	log:
		os.path.join(os.getcwd(), "results/logs/pool_chromosomes/{nchrom}.log")
	shell: 
		'''
		cp {input[1]} {input[0]}/adaptive_gls.csv
		cd {params.genotypooler}/runtools
		python3 parallel_pooling.py {params.cwd}/results/data/{wildcards.nchrom} STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.vcf.gz {input[0]} STU.Chr{wildcards.nchrom}.SNPs.pruned.sorted.pooled.vcf.gz {params.nb_cores} > {log}
		cd {params.cwd}
		rm -r {input[0]}
		'''

