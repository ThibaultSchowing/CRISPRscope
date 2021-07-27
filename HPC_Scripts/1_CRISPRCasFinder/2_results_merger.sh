#!/bin/bash

#SBATCH -o slurm_output/merger-%j-output.txt
#SBATCH -e slurm_output/merger-%j-error.txt

#SBATCH --job-name="merger"
#SBATCH --time=1-00:00:00
#SBATCH --mem=25G
#SBATCH --partition=pall
#SBATCH --export=NONE
#SBATCH --ntasks=1

#======================================================================
#======================================================================
#
# 2_results_merger.sh
#
# For every organism in the defined Parsed_Output/DATE/SOURCE_TYPE
# 	- Create a R tibble for spacers and repeats
# 
# For all organisms
# 	- Merge the created RDS and filter (evidence lvl > 3)
#	- Export Spacers and Repeats for clustering
#	- Clustering: cd-hit / parsing / merging to main tibble
#	- Cas subtype inference: cctyper / parsing / merging
# Inputs variables
#	- CRISPRscope.conf : only for DATE and SOURCE_TYPE
# 			     ./Parsed_output is fixed here (?)
#
# Called scripts: 
#	- 3_results_merger_1.R  
# 	- 4_results_merger_2.R 
#	- 5_clustering.sh
#	- 6_cluster_parser.py
#	- 7_cluster_merger.R  
# 	- 8_cctyper_merger.R
#
#
#======================================================================
#======================================================================

module load vital-it/7
module rm R/3.6.1;
module add R/3.5.1;

mkdir -p slurm_output_merger
rm slurm_output_merger/*

# NEXT STEPS

#=====================================================
# R script V1 -> For each organism Get main dataset
# Get the repeats only after quality filtering (easier from tibble)
#=====================================================

# 

#=====================================================
# Using conda get Cas subtype with Repeattyper
#=====================================================


#=====================================================
# Using the repeat and spacers output files, cluster and import
#=====================================================




# Get the paths, can be executed in any sbatch to have access to the paths and variables
source CRISPRscope.conf



echo "[DEBUG] 2_results_merger.sh"
#[TEST]
#TODO: remove the TEST overwrite before using on real data
#DATE=TEST
#SOURCE_TYPE=TEST


# Parsed results are in : 
# /Parsed_Output/${DATE}/${SOURCE_TYPE}/${ORGANISM}

OutputRmerger=./R_Output
ParsedOutput=./Parsed_Output/${DATE}/${SOURCE_TYPE}

echo "[DEBUG] OutputRmerger: ${OutputRmerger}"
echo "[DEBUG] ParsedOutput: ${ParsedOutput}"


mkdir -p ${OutputRmerger}
rm ${OutputRmerger}/*

for ORGANISM in $(ls ./Parsed_Output/${DATE}/${SOURCE_TYPE}/)
do
	echo ${ORGANISM}
	# Pass variable to R script: https://www.r-bloggers.com/2015/02/bashr-howto-pass-parameters-from-bash-script-to-r/
	
	spacer_file=./Parsed_Output/${DATE}/${SOURCE_TYPE}/${ORGANISM}/${ORGANISM}_PerSpacer_CRISPR_v1.csv
	repeat_file=./Parsed_Output/${DATE}/${SOURCE_TYPE}/${ORGANISM}/${ORGANISM}_DR_v1.fasta
	organism_folder=./Parsed_Output/${DATE}/${SOURCE_TYPE}/${ORGANISM}
	
	# Read the parsed output files of CRISPRCasFinder into tibbles and save the individual files
	# for each organism. (merged in the next loop)
	# UNCOMMENT BEFORE RUNNING
	Rscript 3_results_merger_1.R ${spacer_file} ${repeat_file} ${SOURCE_TYPE} ${organism_folder} ${ORGANISM}
	
		
done



# For each organism (same loop)
# 	read the rds and merge it into a proper big ass dataframe (one for spacer one for repeat) and store it in the output folder ${OutputRmerger}

# Give paths of: Parsed output (for loop)
#		 OutputRmerger (store output dataframes)
#	 	 spacer and repeat files that will be used for clustering and CRISPRCasTyper
# 

filename_output_spacers="CRISPRscope_spacers.fasta"
filename_output_repeats="CRISPRscope_repeats.fasta"

#-----------------------------
# Merging and quality filtering
# Results are in ./R_Output/
#-----------------------------

# UNCOMMENT BEFORE RUNNING
Rscript 4_results_merger_2.R ${ParsedOutput} ${OutputRmerger} ${filename_output_spacers} ${filename_output_repeats} ${DATE} ${SOURCE_TYPE} # loop inside rscript 


# At this point:
# in OutputRmerger: Main dataframe with spacers and repeats
#		    fasta of spacers and repeats
#
# The main dataframe is in rds format to improve loading time with R.

CRISPRscope_main=${OutputRmerger}/${DATE}_${SOURCE_TYPE}_CRISPRscope_tibble.rds



#======================================================================
#======================================================================
# CLUSTERING WITH CD HIT
# 
# The bash script launches cd-hit 
# The python script parses the results and create a csv (Cluster_xxxx)
# 
# The two csv (spacers and repeats) are then joined to the main dataset
#   with R
#
#======================================================================
#======================================================================

spacers_fasta=${OutputRmerger}/${DATE}_${SOURCE_TYPE}_${filename_output_spacers}
repeats_fasta=${OutputRmerger}/${DATE}_${SOURCE_TYPE}_${filename_output_repeats}

if [ -f "${spacers_fasta}" ] && [ -f "${repeats_fasta}" ]; then
	echo "[TEST-OK] ${spacers_fasta} exists " 
	echo "[TEST-OK] ${repeats_fasta} exists"
	
	
	# Clustering Spacers
	./5_clustering.sh ${spacers_fasta} ${OutputRmerger}
	
	spacers_cluster_csv=${OutputRmerger}/Clusters_${DATE}_${SOURCE_TYPE}_${filename_output_spacers%.fasta}.csv
	
	
	# Clustering repeats
	./5_clustering.sh ${repeats_fasta} ${OutputRmerger}
		
	repeats_cluster_csv=${OutputRmerger}/Clusters_${DATE}_${SOURCE_TYPE}_${filename_output_repeats%.fasta}.csv
	
	
	# In R script: load the main dataframe and join the cluster data
	
	echo "[DEBUG] Launch R cluster merger with command: Rscript 7_cluster_merger.R ${CRISPRscope_main} ${spacers_cluster_csv} ${repeats_cluster_csv}"
	
	Rscript 7_cluster_merger.R ${CRISPRscope_main} ${spacers_cluster_csv} ${repeats_cluster_csv}
	

	
	
else
	echo "[TEST-FAILED] ${spacers_fasta} does not exist"
	echo "[TEST-FAILED] ${repeats_fasta} does not exist"
	echo "Skip clustering. "
fi


#======================================================================
#======================================================================
#
# Cas-Subtype inference with CRISPRCasTyper (RepeatTyper)
#
# WARNING: as conda always asks for input when managing envs
#          the cctyper environment should already be installed
#
#
#======================================================================
#======================================================================


# Install the conda environment
#conda env remove --name cctyper
#conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper

echo "============================"
echo "CCTYPER"
echo "============================"

echo "[INFO]: Number of unique repeats: $(cat ${OutputRmerger}/${DATE}_${SOURCE_TYPE}_unique_repeats.tsv | wc -l)"

mkdir -p CCTYPER_OUTPUT_TMP

echo "[DEBUG]: source bashrc and activate env. cctyper"
source ~/.bashrc
conda activate cctyper > ./CCTYPER_OUTPUT_TMP/conda_activate_warning.txt
repeatType ${OutputRmerger}/${DATE}_${SOURCE_TYPE}_unique_repeats.tsv > ./CCTYPER_OUTPUT_TMP/cctyper_raw_output.tsv
conda deactivate

echo "[DEBUG]: Repeat Typer done"
echo "[DEBUG]: Start parsing and merging CCTYPER output. "

# Parsing the CCtyper output
cctyper_output_tsv=./CCTYPER_OUTPUT_TMP/cctyper_raw_output.tsv
Rscript 8_cctyper_merger.R ${CRISPRscope_main} ${cctyper_output_tsv}

rm -r CCTYPER_OUTPUT_TMP

echo "[INFO]: Now ${CRISPRscope_main} contains an up to date dataset for all the genomes in DATE/SOURCE_TYPE"
echo "[END]" 



















