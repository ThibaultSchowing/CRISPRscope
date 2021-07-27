#!/bin/bash

#SBATCH -o slurm_output/CCFDB_align-out.txt
#SBATCH -e slurm_output/CCFDB_align-err.txt

#SBATCH --job-name="AL_CCFDB"
#SBATCH --mem=20G
#SBATCH --partition=pshort


#=========================================================================================
# Compare the Direct Repeats to the CRISPRCasFinder Database (CCFDB)
# and extract the Cas type and species
# 
#  
#=========================================================================================

module add vital-it/7
module load UHTS/Analysis/crass/1.0.1;
module load UHTS/Analysis/sratoolkit/2.10.7;

source CRISPRscope_meta.conf
source ~/miniconda3/etc/profile.d/conda.sh

BIOPROJECT=$1

echo "${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr"

conda activate pycrispr

python3 4a_Align_CCFDB.py ${BIOPROJECT} ${CCFDB_match} ${CCFDB_mismatch} ${CCFDB_open} ${CCFDB_extend} ${CCFDB_threshold}

echo "Now with the CRISPRmap dataset"

python3 5a_Align_CRISPRmap.py ${BIOPROJECT} ${CCFDB_match} ${CCFDB_mismatch} ${CCFDB_open} ${CCFDB_extend} ${CCFDB_threshold}



conda deactivate

















