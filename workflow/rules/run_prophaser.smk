import os
import subprocess
import yaml
from snakemake.utils import min_version
min_version("5.3.0")


# =========================================================================================================
#     Run prophaser
# =========================================================================================================
		
rule run_prophaser:
	"""
	Imputation with prophaser for chromosome 1A only.
	Create output directory for ending the workflow without waiting for the sbatch jobs to start.
	"""
	input:
		"opt/prophaser/main.o"
	output:
		directory("results/data/1/prophaser")
	params:
		studypop_size = 496  # TODO: config parameter
	shell:
		"""
		mkdir results/logs/run_prophaser
		sbatch --array=1-{params.studypop_size} --output results/logs/run_prophaser/slurm-%A_%a.out workflow/scripts/run_workflow_prophaser.sbatch
		mkdir results/data/1/prophaser/
		"""

	
