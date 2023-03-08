#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p node
#SBATCH -C mem128GB
#SBATCH -n 16
#SBATCH -t 00:07:00
#SBATCH -J prophaser

# Location: poolimputeSNPs/
# Run the workflow
apptainer exec container.sif micromamba run -n base bash workflow/scripts/run_workflow_prophaser.sh
