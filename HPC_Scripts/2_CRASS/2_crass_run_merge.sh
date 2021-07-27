#!/bin/bash

#SBATCH -o slurm_output/crass-output-%j.txt
#SBATCH -e slurm_output/crass-error-%j.txt

#SBATCH --job-name="CRASS"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=20G
#SBATCH --partition=pall

#=========================================================================================
#
#   Executes CRASS on the downloaded metagenome (by Bioproject number)
#   Merges the result files into one .crispr file per bioproject
#
#
#=========================================================================================


module add vital-it/7
module load UHTS/Analysis/crass/1.0.1;



module load UHTS/Analysis/sratoolkit/2.10.7;

BIOPROJECT=$1

source ./CRISPRscope_meta.conf

# Directory containing the BioProjects folders

RES=${RESDIR}

echo "_____________CRASS____________"
echo " "
echo "Result directory: ${RESDIR} "
echo "Bioproject: ${BIOPROJECT}"


# Check if Data exist
if [[ -d "${DATA}/${BIOPROJECT}" ]]; then
	echo "Data available at ${DATA}/${BIOPROJECT} -> continuing"
else
	echo "Data not found at ${DATA}/${BIOPROJECT} -> quitting."
	echo "! ABORT: no data !"
	exit 0
fi


DATA_PRJ=$(ls $DATA/$BIOPROJECT)

#TODO: NAME OF FASTQ FILE NOT EQUAL TO DIRNAME !!! CHECK FOR CUSTOM VINCENT DATA -> DO SOMETHING FOR ALL PROJECTS
# arr = list of sequence names
echo "List of entries in ${BIOPROJECT}"
echo "-> ${DATA_PRJ}"
read -a arr <<< $DATA_PRJ

echo "================================================================="
echo "Entries available for the selected Bioproject or directory:"
echo "================================================================="
echo " ${#arr[@]}"

echo "======="

###shopt -s nullglob
###array=($DATA$BIOPROJECT/*/)
###shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later
###echo "${array[@]}"  # Note double-quotes to avoid extra parsing of funny characters in filenames



if (( ${#arr[@]} == 0 )); then
	echo "ERROR: Array length is 0. No sequence folder found" >&2
	echo "ERROR: Array length is 0. No sequence folder found"
fi


###test=(ls -d -- $DATA$BIOPROJECT/*/)

# For all the file (DIRECTORIES!!) of the BioProject data, execute crass and save the results (merge .crispr files later)
for f in ${arr[@]} 
do
	#FILES=${DATA}/${BIOPROJECT}/${f}/${f}*.fastq.gz
	
	#for fastfile in $FILES
	#do
	#	OUT=${RES}/${BIOPROJECT}/${f}
	#	mkdir -p ${OUT}
	#	echo "Output directory: ${OUT}"
	#	echo "Input file: ${fastfile}"
	#	if [[ -f "${fastfile}" ]]; then
	#		echo "File exists"
	#	else
	#		echo "File ${fastfile} not found"
	#	fi
	#	
	#	# If CRASS has already been executed, skip.
	#	if [[ -f "${OUT}/crass.crispr" ]]; then
	#		echo "[WARNING] file ${OUT}/crass.crispr already exists for ${f}.. skipping."
	#		continue
	#	else
	#		
	#		echo "Start CRASS for $f"
	#		crass -o ${OUT} ${fastfile}
	#		echo "CRASS done for ${f}"
	#	fi
	#done
	
	OUT=${RES}/${BIOPROJECT}/${f}
	rm -r ${OUT}
	mkdir -p ${OUT}
	echo "========================"
	echo "Output directory: ${OUT}"
	#FQ_FILES="${DATA}/${BIOPROJECT}/${f}/*.fastq.gz"
	#echo ${FQ_FILES}
	
	# For all the fastq files in each read directories execute crass
	# Use globbing -> the variable contains the entire path to the file
	#for fastfile in $(ls ${DATA}/${BIOPROJECT}/${f}/*.fastq.gz)
	for fastfile in ${DATA}/${BIOPROJECT}/${f}/*.fastq.gz
	do
		[ -e "$fastfile" ] || continue

		echo ""	
		echo "Input fastq.gz: -> ${fastfile}"
		echo "=Launching crass="
		crass -o ${OUT} ${fastfile}
			
	done

done





# List all the crass.crispr files in the Bioproject subdirectories (SRR or SRA entries)
#crispr_lst=$(find ${RES}/${BIOPROJECT} -name 'crass.crispr')

#===========================================================
#  MERGING = MOVED IN PARSING STEP
#===========================================================

#echo "Merging the following .crispr files: "
#
#count=0
#for crsp in ${crispr_lst[@]}
#do
#	
#	echo $crsp
#	count=$(($count + 1))
#	
#	#
#	# Extract the fasta files independently for each entry
#	#
#	
#	
#	
#	
#	
#done
#
#echo "Destination: ${RES}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr"
#
## Test if the number of .crispr files is bigger than 1
#
#if (( ${count} > 1 )); then
#	echo "More than one .crispr file... merging."
#	crisprtools merge -s -o ${RES}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr $crispr_lst
#else
#	echo "Only one .crispr file... moving ( number of elements: $count x${#crispr_lst[@]} )"
#	echo "file: ...\n${crispr_lst}"
#	echo "..."
#	cp ${crispr_lst} ${RES}/${BIOPROJECT}
#	mv ${RES}/${BIOPROJECT}/${crispr_lst} ${RES}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr
#fi
#
#

#echo Crass ran and summary files merged for Bioproject $BIOPROJECT



























