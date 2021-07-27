#!/bin/bash

#SBATCH -o slurm_output/preprocess-output.txt
#SBATCH -e slurm_output/preprocess-error.txt

#SBATCH --job-name="preproc"
#SBATCH --time=3:00:00
#SBATCH --mem=5G
#SBATCH --partition=pshort
#SBATCH --export=NONE




# Load modules

module load vital-it/7

#
# Preprocess the files before running the main script
# Files format / naming or origin depends on their sources
# Project was first based only on Dialact but now also on NCBI assemblies -> we need short and unique names for each assembly
# 
#

ORGANISM=$1

if [ "$#" -ne 1 ]; then
	echo "Illegal number of parameters"
	exit
fi


# SOURCE IN MANUAL
#DATE=20201201
#DATE=20210406
#DATE=20210510
#DATE=20210511
DATE=20210514

SOURCE=/data/projects/p539_crisprscope/TMP_SOURCE_NCBI/${DATE}
#SOURCE=/data/projects/p539_crisprscope/TMP_SOURCE_DIALACT/${DATE}
SOURCE_TYPE=NCBI

# Root of the project. Needed to bind within the Singularity container. Needs to contain the scripts and
# all the data.

ROOTDIR=..

# Contains the assemblies. The same structure as SOURCE will be kept (folders by organisms containing
# the related assembles by strain). The fasta files need to be adapted for CRISPRCasFinder and thus
# they will be copied and modified here.

DATA=../0_data/assemblies



# Input directory:
# Assemblies will be copied here here in specific organism folders (folder created bellow)
INDIR=${DATA}/${SOURCE_TYPE}/${DATE}/${ORGANISM}

# Directory where this script is stored (wow so useful !)
SCRIPTDIR=.


# Here we only need the INDIR file
# Copied assemblies, will be modified (format)
mkdir -p ${INDIR}


# Create the directories
#source CRISPRscope_init.sh $ORGANISM

echo "Configuration: "
echo ""
echo "Source type: ${SOURCE_TYPE}"
echo "Source: ${SOURCE}"
echo "Date: ${DATE}"
echo "Indir: ${INDIR}"




if [[ ${SOURCE_TYPE} == "NCBI" ]]; then
	
	
	FILENAME_PATTERN_NCBI_after="^[0-9]+.[0-9].fna$"
	FILENAME_PATTERN_NCBI="^GCF_[0-9]+_*.fna$"

        # Maybe for NCBI, just check the two first fields and the .fna in the end.
	
	# Data from NCBI: check the file names and rename with shorter/unique ID
	# Take the second field of the filename (format: GCF_008690745.1_ASM869074v1_genomic.fna)
	
	echo "Source type: ${SOURCE_TYPE}"
	
	# List files in TMP_SOURCE_NCBI
	
	NBFILES=$(ls -1q "${SOURCE}/${ORGANISM}"/ | wc -l)
	echo "There are ${NBFILES} files in the source directory."
	if [ "${NBFILES}" -eq "0" ]; then
		echo "There are no files to copy. Verify your downloads."
		echo "Exit."
		exit 0
	fi
	
	# For each file:
	# Check if it already exists in ${INDIR} 
	
	# copy and shorten the rename
	
	
	for file in "${SOURCE}/${ORGANISM}"/*.fna; do
		echo "----------------------"
		echo "Checking file: ${file}"
		
		#TODO: here we don't check the pattern of the filename. Should be verified but we can count on the consistency of NCBI.
		# IF not like pattern (starts with ^GCF_[0-9]+_*.fna
		if [[ ! ${file##*/} =~ $FILE_NAME_PATTERN_NCBI ]]; then
			echo "NCBI filename does not match the subseting pattern."
			echo "Skip this file."
			continue
		fi
		
		
		
		# Take shorter filename (substring, field 2, sep _
		filename=${file##*/}
		echo "Filename: ${filename}"
		
		filerenamed=GCF.$(echo ${filename} | cut -d'_' -f 2)
		filerenamed+=".fna"
		echo "Filerenamed: ${filerenamed}"
		
		
		
		
		# IF files has already been preprocessed, skip.
		echo "Check if ${INDIR}/${filerenamed} exists."
		
		if [[ -f "${INDIR}/${filerenamed}" ]]; then
			# file exists, do something
			echo "${filerenamed} already exists in ${INDIR}"
			echo "===> ${INDIR}/${filerenamed##*/}"
		else
			echo "File doesn't exist in destination directory.. copying."
			echo "cp ${file} ${INDIR}/${filerenamed}"
			cp ${file} ${INDIR}/${filerenamed}
			echo "Done for file ${filerename}"
		fi
	
	
	
	done
	
	
	
elif [[ ${SOURCE_TYPE} == "DIALACT" ]]; then
	
	FILENAME_PATTERN_DIALACT="^[a-zA-z][0-9]+[-][ip][0-9][-][0-9].fna$"
	
	# Data from dialact: check the file names (mean recent PGAP annotated etc) 
	
	echo "Source type: ${SOURCE_TYPE}"

        # List files in TMP_SOURCE_DIALACT

        NBFILES=$(ls -1q "${SOURCE}/${ORGANISM}/"*.fna | wc -l)
        echo "There are ${NBFILES} files in the source directory."
        if [ "${NBFILES}" -eq "0" ]; then
                echo "There are no files to copy. Verify your downloads."
                echo "Exit."
                exit 0
        fi

	
	for file in "${SOURCE}/${ORGANISM}"/*.fna; do
		echo "----------------------"
		echo "Checking file: ${file}"
	
		filename=${file##*/}
                echo "Filename: ${filename}"
		
		
		
		# IF not like pattern (starts with ^GCF_[0-9]+_*.fna
                if [[ ! ${file##*/} =~ $FILE_NAME_PATTERN_DIALACT ]]; then
                        echo "DIALACT filename does not match the subseting pattern."
                        echo "Skip this file"
			continue
                fi
		
		
		
		# IF files has already been preprocessed, skip.
		echo "Check if ${INDIR}/${filename} exists."
		if [[ -f "${INDIR}/${filename}" ]]; then
			# file exists, do something
			echo "${filename} already exists in ${INDIR}"
			echo "===> ${INDIR}/${filename##*/}"
		else
			echo "File doesn't exist in destination directory.. copying."
			echo "cp ${file} ${INDIR}/${filename}"
			cp ${file} ${INDIR}/${filename}
			echo "Done for file ${filename}"
		fi
		
		
		
		
		
		
		
		
		
	done
	
	
	
	
	
# Unexisting source / add case here if you add a new source. 	
else
	echo "Unknown data source/type"
	exit 0
fi



echo "Pre-processing done. END."
exit 1















# Copy files from their original location (buffer) to the input directory
# Do not copy if files are already in the assemblies directory
# Do not copy files not matching the regex.


for file in "${SOURCE}/${ORGANISM}"/*.fna; do

        #echo "${INDIR}/${file##*/}"

        if [[ -f "${INDIR}/${file##*/}" ]]; then
                # file exists, do something
                echo "${file} already exists"
                echo "===> ${INDIR}/${file##*/}"
                echo " "
        else
                # Check if file matches regex
                if [[ ${file##*/} =~ $FILE_NAME_PATTERN ]]
                then
                        echo "${file}"
                        echo "${file##*/}"
                        echo "File OK .. copying ${file} in ${INDIR}"
                        # Copy the file from the main place to the 0_data/assemblies/organism folder.
                        cp ${file} $INDIR
                else
                        echo "/!\ File format not matching /!\ "
                        echo "File: ${file##*/} will not be copied!"
                        echo "Quitting program now. \nReason: fna file name not matching last naming convention."
                        exit
                fi
        fi
done

echo "--------------------------"
echo "Copy: Done !"

