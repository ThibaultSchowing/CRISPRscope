#!/bin/bash

#SBATCH -o slurm_output/fetcher-output.txt
#SBATCH -e slurm_output/fetcher-error.txt

#SBATCH --job-name="DataFetch"
#SBATCh --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH --partition=pall



module load vital-it/7

# As prefetch from sra-tools is unreliable, an alternative solution has to be inplemented.
# Issue: https://github.com/ncbi/sra-tools/issues/276

# Note that prefetch stoped working without reason at a random time.

# As alternative, enaDataGet will be used. Info at: https://github.com/enasequence/enaBrowserTools

# Don't forget to do:

# alias enaDataGet=/data/projects/p539_crisprscope/enaBrowserTools/python3/enaDataGet
# alias enaGroupGet=/data/projects/p539_crisprscope/enaBrowserTools/python3/enaGroupGet


# Usage: enaDataGet --help
# for bioproject use: enaGroupGet --help


BIOPROJECT=$1

source CRISPRscope_meta.conf

force_new=false


#=============================================================
# Verifications
#=============================================================

#Note: doesn't check for each SRR entry but for the whole bioproject (empty or not) 
#      be warned, it takes time. 

#
# Creates the data directory for the bioproject if needed or skip
#
if [[ -d "${DATA}/${BIOPROJECT}" && "$force_new" = true  ]] ; then
        echo "-----"
        echo "Data directory already exists..."
	echo "Forcing new download..."
        echo "Deleting all entries..."
        rm -rf ${DATA}/${BIOPROJECT}
        echo "Done!"
        echo "Creates output directory ${DATA}/${BIOPROJECT}"
        mkdir -p ${DATA}/${BIOPROJECT}
        echo "Done!"
        echo "-----"
elif [[ -d "${DATA}/${BIOPROJECT}" && "$force_new" = false ]] ; then
        echo "-----"
        echo "Data directory already exists..."
        echo "Keeping old data..."
        echo ""
	
	if [[ "$(ls -A ${DATA}/${BIOPROJECT})" ]]; then
		# directory is not empty 
	        echo "Data already present, resume download (hope it's automatic)"
		#echo "-----> exit"
		#exit 1
	else
		echo "Directory ${DATA}/${BIOPROJECT} exists but is empty -> download data"
		rm -r ${DATA}/${BIOPROJECT}/*
	fi

elif [[ ! -d "${DATA}/${BIOPROJECT}" ]] ; then
        echo "-----"
	echo "Data directory doesn't exist..."
        echo "Creating output directory..."
        mkdir -p ${DATA}/${BIOPROJECT}
        echo "Done!"
        echo "-----"
fi






#=============================================================
# Download the files to the DESTINATION directory
#=============================================================


# without -d -> ./


# Let's go for the long download
echo "Start download for Bioproject ${BIOPROJECT} to the folder: ${DATA}..."
STARTTIME=$(date +%s)

${ENA}/python3/enaGroupGet --format fastq --dest ${DATA} $BIOPROJECT

ENDTIME=$(date +%s)


echo "Done in $(($ENDTIME - $STARTTIME)) seconds. Have a nice day."












