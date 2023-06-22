#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 00:07:00
#SBATCH -J poolimputeSNPs

# Location: poolimputeSNPs/
# Run only imputation with prophaser
sbatch --array=1-496 --output results/logs/run_prophaser/slurm-%A_%a.out workflow/scripts/run_workflow_prophaser.sbatch
