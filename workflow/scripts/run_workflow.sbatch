#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:10:00
#SBATCH -J poolimputeSNPs

# Location: poolimputeSNPs/
# Run the workflow
bash workflow/scripts/run_workflow.sh
