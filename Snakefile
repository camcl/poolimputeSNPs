"""
Main Snakemake file for making Snakemake to chain the intermediate snakefiles.
Output driven: last output of the last intermediate Snakemake file to be run. Snakemake will figure out how to chain for making the pipeline to work.
"""

# =================================================================================================

from snakemake.utils import min_version

# Make Sure Minimun Snakemake version
min_version("5.3.0")

# =========================================================================================================
#     Setup Config and Report
# =========================================================================================================

# Load config file
configfile: "config/config.yml"

# validate(config, schema="../schemas/config.schema.yaml")

# Description of the workflow can be found in the final report
#report: "report/workflow.rst"

# container: "docker://continuumio/miniconda3"

# =========================================================================================================
#     Load Rules
# =========================================================================================================

include: "workflow/rules/install_genotypooler.smk"
include: "workflow/rules/prepare_stu_data_55K.smk"     
include: "workflow/rules/prepare_pnl_data.smk"
include: "workflow/rules/eda.smk"
include: "workflow/rules/pool_chromosomes.smk"
include: "workflow/rules/genetic_maps.smk"

# =========================================================================================================
#     The `onstart` Checker
# =========================================================================================================

'''
onstart:
    try:
        print("Checking if all required files are provided...")
        important_files = [ config["samples"],
                            config["phenotypes"] ]
        for filename in important_files:
            if not os.path.exists(filename):
                raise FileNotFoundError(filename)
    except FileNotFoundError as e:
        print("This file is not available or accessible: %s" % e)
        sys.exit(1)
    else:
        print("\tAll required files are present!")
'''


# =========================================================================================================
#     Target Outputs
# =========================================================================================================

rule all:
    input:
        "results/rulegraph.png",
        "results/jobgraph.png"
        

rule generate_workflow_graphs:
	"""
	Generate a rulegraph for the entire workflow.
	"""
	input:
		# expand("results/data/{nchrom}/STU.Chr{nchrom}.SNPs.pruned.sorted.pooled.{ext}", nchrom=list(config["chromosomes"]["prefix"].values()), ext=["vcf.gz", "vcf.gz.csi"])
		expand("results/data/{nchrom}/{nchrom}_interpolated_wheat_map", nchrom=list(config["chromosomes"]["prefix"].values())) # output from genetic_maps.smk
	output:
		rulegraph = "results/rulegraph.png",
		jobgraph = "results/jobgraph.png"
	shell:
		"""
		snakemake --rulegraph | dot -Tpng > {output.rulegraph}
		snakemake --dag | dot -Tpng > {output.jobgraph}
		"""

# =========================================================================================================
#     Success and Failure Messages
# =========================================================================================================

onsuccess:
	print("Success! The Snakemake workflow is completed.")

onerror:
	print("Error! The Snakemake workflow aborted.")
