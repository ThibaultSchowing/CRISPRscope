#!/bin/bash

#SBATCH -o slurm_output/master-output-%j.txt
#SBATCH -e slurm_output/master-error-%j.txt
#SBATCH --job-name="crscp_meta"
#SBATCH --ntasks=4

#SBATCH --time=10-24:00:00
#SBATCH --mem=20G
#SBATCH --partition=pall






BIOPROJECT=$1

source CRISPRscope_meta.conf

source ~/miniconda3/etc/profile.d/conda.sh

echo "CRISPRscope Meta - Master script"
echo "--------------------------------"
echo " "
echo "Bioproject: ${BIOPROJECT}"




# CCF1=$(sbatch --wait 1_CRISPRCasFinder_task_launcher.sh ${ORGANISM})
# echo ${CCF1}
# CCF1_JID=${CCF1##* }




#====================================================================
# Metagenome fetcher
#
# -> if the data doesn't come from BioProject check conf file (NOFETCH)
#====================================================================

echo "Launching the Bioproject data fetcher job... "
if [[ "${NOFETCH}" == false ]] ; then
	FETCHER=$(sbatch --wait 1_bioproject_metagenome_fetcher.sh ${BIOPROJECT})
	FETCHER_JID=${FETCHER##* }
	echo "Job launched -> FETCHER_JID:${FETCHER_JID}"


#====================================================================
# CRASS - assemble and merge results
#====================================================================

	echo "Launching CRASS on the metagenomic data or waiting for dependency on job ${FETCHER_JID}"
	CRASSJOB=$(sbatch --wait --dependency=afterany:${FETCHER_JID} 2_crass_run_merge.sh ${BIOPROJECT})
	CRASSJOB_JID=${CRASSJOB##* }

	echo "Job launched -> CRASSJOB_JID:${CRASSJOB_JID}"
else
	# If the data doesn't need to be fetched (not SRA related), 
	
	echo "Launching CRASS on the metagenomic data without depenency"
	CRASSJOB=$(sbatch --wait 2_crass_run_merge.sh ${BIOPROJECT})
        CRASSJOB_JID=${CRASSJOB##* }

fi





#====================================================================
# Extract data from the .crispr file with sratoolkit
#====================================================================


echo "Launching sequence and stats extraction with SRAtoolkit or waiting on job ${CRASSJOB_JID}"
ASDF=$(sbatch --wait --dependency=afterany:${CRASSJOB_JID} 3_parsing_crass.sh ${BIOPROJECT})
ASDF_JID=${ASDF##* }

echo "Job launched -> PARSING_JID:${ASDF_JID}"



#====================================================================
# Align the repeats with the repeats databases CCFDB and CRISPRmap
#====================================================================

echo "Launching sequence matching with CRISPRCasFinder Database and CRISPRmap database. "


QWERT=$(sbatch --wait --dependency=afterany:${ASDF_JID} 4_get_castype_CCFDB.sh ${BIOPROJECT})
QWERT_JID=${QWERT##* } 











echo "End of CRISPRscope Meta Job launcher"











