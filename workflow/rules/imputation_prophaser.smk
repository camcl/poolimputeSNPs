import subprocess
import yaml
from snakemake.utils import min_version
min_version("5.3.0")


# =========================================================================================================
#     Clone repository from GitHub, checkout branch, retrieve commit, and fill in the config file with these infos
# =========================================================================================================
	
rule clone_compile_prophaser: 
	"""Clone the repository, write commit reference to the configuration file and compile the code."""
	input:
		expand("results/data/{nchrom}/{nchrom}_interpolated_wheat_map", nchrom=list(config["chromosomes"]["prefix"].values()))
	output:
		"opt/prophaser/main.o"
	params:
		url="https://github.com/scicompuu/prophaser.git"
	log:
		"results/logs/imputation_prophaser/clone_compile_prophaser.log"
	shell:
		"""
		sed -i '/commit_hash_prophaser:/d' config/config.yml
		cd opt/
		git clone https://github.com/scicompuu/prophaser.git
		cd prophaser
		git checkout multilevel
		rm makefile
		cp ../../makefile makefile
		make
		hash_commit=$(git rev-parse --short HEAD)
		echo "commit hash prophaser: " > ../../{log}
		echo "  $hash_commit" >> ../../{log}
		echo "" >> ../../config/config.yml
		echo "commit_hash_prophaser: $hash_commit" >> ../../config/config.yml
		cd ../..
		"""
		
	
