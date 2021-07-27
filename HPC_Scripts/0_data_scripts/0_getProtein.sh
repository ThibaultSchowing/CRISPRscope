#!/bin/bash

#SBATCH -o getgen-out.txt
#SBATCH -e getgen-err.txt

#SBATCH --job-name="GetGen"
#SBATCH --mem=20G
#SBATCH --partition=pshort
#SBATCH --time 3:00:00


module add vital-it/7
module add Blast/ncbi-blast/2.10.1+;




##-----------------
##NCBI
##-----------------

TARGET_DIR=/data/projects/p539_crisprscope/TMP_SOURCE_NCBI

# /data/projects/p539_crisprscope/0_data/assemblies/NCBI/
date=20210514
for species in $(cat ./SpeciesNames.txt)
#for species in $(cat $mainPath/02_script/names.txt)
do

	echo ==========================
	echo ${species}
	echo ==========================


	#species=Streptococcus_thermophilus


	##log File

	mkdir -p ${TARGET_DIR}/${date}/logs

	speciesSpace=$(echo $species | awk -F "_" 'BEGIN{OFS=" "} {print $1,$2}')

	curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' > ${TARGET_DIR}/${date}/logs/logfile_ncbi_genomes_all.txt
	curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' | gawk 'BEGIN{FS="\t";OFS="\t"} /^#/ {next} {print $0}' | \
	grep "${speciesSpace}" > ${TARGET_DIR}/${date}/logs/logfile_ncbi_genomes_${species}.txt

	cut -f 20 ${TARGET_DIR}/${date}/logs/logfile_ncbi_genomes_${species}.txt | \
	cut -d '/' -f 10 > ${TARGET_DIR}/${date}/logs/genome_names_${species}.txt
	curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt' > ${TARGET_DIR}/${date}/logs/README_assembly_summary.txt

	##==============
	##IMPORTANT!!! 
	##CHANGE PROTEIN/GENOMIC AND FILE FORMAT
	##==============

	##download
	curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' | \
	gawk 'BEGIN{FS="\t";OFS="\t"} /^#/ {next} {print $8 , $20}' | \
	grep "${speciesSpace}"  | cut -f 2 | \
	sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)|\1\2\/\2_genomic.fna.gz|' > ${TARGET_DIR}/${date}/logs/downloadFile_${species}.txt


	mkdir -p ${TARGET_DIR}/${date}/${species}/
	cd ${TARGET_DIR}/${date}/${species}/
	
	# Take only the first 500 genomes
	head -n 500 ${TARGET_DIR}/${date}/logs/downloadFile_${species}.txt > ${TARGET_DIR}/${date}/logs/downloadFile_short_${species}.txt
	rm ${TARGET_DIR}/${date}/logs/downloadFile_${species}.txt
	wget -i ${TARGET_DIR}/${date}/logs/downloadFile_short_${species}.txt 

	gunzip *.fna.gz
	#rm *.fna.gz*

done




echo "ALL GENOMES AVAILABLE DOWNLOADED"



