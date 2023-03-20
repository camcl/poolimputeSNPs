#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p node
#SBATCH -C mem128GB
#SBATCH -n 16
#SBATCH -t 00:30:00
#SBATCH -J poolimputeSNPs

# Location: poolimputeSNPs/
# Run the workflow
bash workflow/scripts/run_workflow.sh
