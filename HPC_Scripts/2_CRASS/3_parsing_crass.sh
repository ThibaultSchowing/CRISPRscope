#!/bin/bash

#SBATCH -o slurm_output/parsing-out-%j.txt
#SBATCH -e slurm_output/parsing-err-%j.txt
#SBATCH --time=90:00:00
#SBATCH --job-name="sratools"
#SBATCH --mem=10G
#SBATCH --partition=pall

#=========================================================================================

module add vital-it/7
module load UHTS/Analysis/crass/1.0.1;
module load UHTS/Analysis/sratoolkit/2.10.7;

source CRISPRscope_meta.conf


BIOPROJECT=$1


#===========================================================
# For each SRR entry, extract fasta + keep metadata
#===========================================================


for D in ${RESDIR}/${BIOPROJECT}/*; do
#    if [ -d "${D}" ] && [[ -f "${D}/crass.crispr" ]]; then
    if [[ -f "${D}/crass.crispr" ]]; then
	echo "Apparently file ${D}/crass.cripr exists"
	# SRR entry without path: ERR12345678
        echo "Extracting information from sample ${D##*/} ... "
	ENTRY=${D##*/}
	# Extract everything from the crass.crispr
	
	
	crisprtools stat -H ${RESDIR}/${BIOPROJECT}/${ENTRY}/crass.crispr > ${RESDIR}/${BIOPROJECT}/${ENTRY}/$BIOPROJECT.${ENTRY}.stats.csv
	
	#=======================================================
	#    Get Direct Repeats out of the merged .crispr file
	#=======================================================

	crisprtools extract -d ${RESDIR}/${BIOPROJECT}/${ENTRY}/crass.crispr > ${RESDIR}/${BIOPROJECT}/${ENTRY}/$BIOPROJECT.${ENTRY}.DR.fasta
	
	
	#=======================================================
	#       Get spacers out of the merged .crispr file
	#=======================================================
	
	crisprtools extract -s ${RESDIR}/${BIOPROJECT}/${ENTRY}/crass.crispr > ${RESDIR}/${BIOPROJECT}/${ENTRY}/$BIOPROJECT.${ENTRY}.Spacers.fasta

	
	#=======================================================
	#    Get flank. seq. out of the merged .crispr file
	#=======================================================


	crisprtools extract -f ${RESDIR}/${BIOPROJECT}/${ENTRY}/crass.crispr > ${RESDIR}/${BIOPROJECT}/${ENTRY}/${BIOPROJECT}.${ENTRY}.FlankingSequences.fasta

	
	
	#
	# Count the number of reads in the related data file 
	#
	
	# TODO: do this in the computation next time
	
	READ_COUNT=$(zcat ${DATA}/${BIOPROJECT}/${ENTRY}/*.fastq.gz | grep "^@" | wc -l) 
	echo "Fasta and stat files extracted from sample ${ENTRY}"
	echo "Read count: ${READ_COUNT}"
	echo "${READ_COUNT}" > ${RESDIR}/${BIOPROJECT}/${ENTRY}/${BIOPROJECT}.${ENTRY}.readcount.txt
	
    fi
done




# List all the crass.crispr files in the Bioproject subdirectories (SRR or SRA entries)
crispr_lst=$(find ${RESDIR}/${BIOPROJECT} -name "crass.crispr")

#===========================================================
#  MERGING 
#===========================================================

echo "Merging the following .crispr files: "
echo "List: ${crispr_lst[@]}"


count=0
for crsp in ${crispr_lst[@]}
do

       echo $crsp
       count=$(($count + 1))


done

echo "Destination: ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr"

# Test if the number of .crispr files is bigger than 1

if [[ ${count} > 1 ]]; then
       echo "More than one .crispr file... merging."
       crisprtools merge -s -o ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr $crispr_lst
elif [[ ${count} == 1 ]]; then
       echo "Only one .crispr file... moving ( number of elements: $count x${#crispr_lst[@]} )"
       echo "file: ...\n${crispr_lst}"
       echo "..."
       cp ${crispr_lst} ${RESDIR}/${BIOPROJECT}
       mv ${RESDIR}/${BIOPROJECT}/${crispr_lst} ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr
else
	echo "Empty list, no crispr file"
fi












#=======================================================
#    Get statistics out of the merged .crispr file
#=======================================================

echo "${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr"

crisprtools stat -H ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr > ${RESDIR}/${BIOPROJECT}/$BIOPROJECT.stats.csv

#=======================================================
#    Get Direct Repeats out of the merged .crispr file
#=======================================================

crisprtools extract -d ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr > ${RESDIR}/${BIOPROJECT}/$BIOPROJECT.DR.fasta


#=======================================================
#       Get spacers out of the merged .crispr file
#=======================================================

crisprtools extract -s ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr > ${RESDIR}/${BIOPROJECT}/$BIOPROJECT.Spacers.fasta


#=======================================================
#    Get flank. seq. out of the merged .crispr file
#=======================================================


crisprtools extract -f ${RESDIR}/${BIOPROJECT}/${BIOPROJECT}_crispr_merged.crispr > ${RESDIR}/${BIOPROJECT}/$BIOPROJECT.FlankingSequences.fasta







