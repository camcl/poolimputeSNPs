#!/bin/bash -l
#SBATCH -A snic2022-22-697
#SBATCH -M snowy
#SBATCH -p node
#SBATCH -C mem128GB
#SBATCH -n 16
#SBATCH -t 03:00:00
#SBATCH -J poolimputeSNPs

sbatch workflow/scripts/beagle_magic_wheat_map.sh
