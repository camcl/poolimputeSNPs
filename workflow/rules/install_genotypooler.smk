import subprocess
import yaml
from snakemake.utils import min_version
min_version("5.3.0")


# =========================================================================================================
#     Clone repository from GitHub, checkout branch, retrieve commit, and fill in the config file with these infos
# =========================================================================================================
	
rule clone_repository: 
	output:
		directory("opt/genotypooler/genotypooler"),
		"opt/genotypooler/genotypooler/runtools/rm_swapped_ref_alt.py"
	params:
		url="https://github.com/camcl/genotypooler.git"
	log:
		"results/logs/install_genotypooler/clone_repository.log"
	shell:
		"""
		sed -i '/commit_hash_genotypooler:/d' config/config.yml
		rm -r opt
		mkdir opt
		cd opt/
		git clone https://github.com/camcl/genotypooler.git
		cd genotypooler 
		git checkout magicwheat-dev
		pip install -r requirements.txt
		hash_commit=$(git rev-parse --short HEAD)
		echo "commit hash genotypooler: " > ../../{log}
		echo "  $hash_commit" >> ../../{log}
		echo "" >> ../../config/config.yml
		echo "commit_hash_genotypooler: $hash_commit" >> ../../config/config.yml
		cd ../..
		"""
	
