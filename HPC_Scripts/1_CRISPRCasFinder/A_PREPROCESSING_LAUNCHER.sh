#!/bin/bash

#SBATCH -o slurm_output/launcher_preprocess-output.txt
#SBATCH -e slurm_output/launcher_preprocess-error.txt

#SBATCH --job-name="lchrprep"
#SBATCH --time=03:00:00
#SBATCH --mem=5G
#SBATCH --partition=pshort
#SBATCH --export=NONE

# Will read the Species names (phase 2 - all species detected in metagenomes)

while IFS="" read -r p || [ -n "$p" ]
do
	#printf '%s\n' "$p"
	#echo "sbatch PREPROCESS_.sh $p"
	sbatch /data/projects/p539_crisprscope/1_ccf_scripts/A_PREPROCESSING_.sh $p
done < PreprocessSpeciesNames.txt




