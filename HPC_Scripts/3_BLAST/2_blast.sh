#!/bin/bash

#SBATCH -o slurm_output/blast-IMG2-output.txt
#SBATCH -e slurm_output/blast-IMG2-error.txt
#SBATCH --time=12-00:00:00
#SBATCH --job-name="BLASTn"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=37
#SBATCH --mem=20G
#SBATCH --partition=pall

# The viral databases have to be prepared with makeblastdb
# -> done in the databases folders


# Goal: nblast the spacers against the viral databases and extract the viral species that matches 
# TODO: also check for PAM on the protospacers. 


module add vital-it/7

module rm Blast/ncbi-blast/2.10.1+;
module add Blast/ncbi-blast/2.9.0+;

module rm Blast/ncbi-blast/2.7.1+;


# makeblastdb -max_file_sz 1GB -in xxxxxxx.fasta -out yyyyyyyy.fasta




#Directory where to copy fasta files of interest
DATA_DIR="./A_data_blast_s"
#DATABASE=../0_data/DB/PhageWeb/PhageWebDB
DATABASE=../0_data/DB/IMGVR/IMG_VR_2020-10-12_5.1/IMGVR

RESULTDIR=./B2_result_blast_IMGVR
RESULTFILE=RES_IMGVR2_results.txt


# List of fasta files of interest to blast:

FASTA_FILES=$(ls ${DATA_DIR})
read -a arr <<< ${FASTA_FILES}

echo "-------------------------------------"
echo "Fasta files for blasting:"
echo "-------------------------------------"
echo ""
echo "${#arr[@]}"

if (( ${#arr[@]} == 0 )); then
    echo "No sequence found" >&2
fi




# For each fasta file blast against the database
for fastafile in ${arr[@]}
do
	echo "blast for file ${fastafile}"
	f=${DATA_DIR}/${fastafile}
	echo "processing file: ${f}"  >&2
	
	blastn -num_threads 37 -max_hsps 1 -max_target_seqs 1 -task blastn-short \
	  -query ${f} \
	  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen sgi sacc sblastnames staxids sscinames" \
	  -db ${DATABASE} \
	  -out ${RESULTDIR}/${fastafile##*/}_IMGVR_nblast_out.txt -evalue 0.01 -word_size 16	
	
	cat ${RESULTDIR}/${fastafile##*/}_IMGVR_nblast_out.txt >> ./${RESULTFILE}
	rm ${RESULTDIR}/${fastafile##*/}_IMGVR_nblast_out.txt
done





