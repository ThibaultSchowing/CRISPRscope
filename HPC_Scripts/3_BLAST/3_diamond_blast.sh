#!/bin/bash

#SBATCH -o slurm_output/diamond-blast-output.txt
#SBATCH -e slurm_output/diamond-blast-error.txt
#SBATCH --time=10-00:00:00
#SBATCH --job-name="BLASTx"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --partition=pall

# The viral databases have to be prepared with makeblastdb
# -> done in the databases folders


# Goal: nblast the spacers against the viral databases and extract the viral species that matches 
# TODO: also check for PAM on the protospacers. 


module add vital-it/7

#module rm Blast/ncbi-blast/2.10.1+;
#module add Blast/ncbi-blast/2.9.0+;
#module rm Blast/ncbi-blast/2.7.1+;
module add SequenceAnalysis/SequenceSearch/diamond/0.9.31;

# makeblastdb -max_file_sz 1GB -in xxxxxxx.fasta -out yyyyyyyy.fasta




#Directory where to copy fasta files of interest
DATA_DIR="./C_data_diamond"
#DATABASE=../0_data/DB/PhageWeb/PhageWebDB
DATABASE=../0_data/DB/IMGVR/IMG_VR_2020-10-12_5.1_Diamond/IMGVR_all_proteins
OUTPUTDIR="D_results_diamond"

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

	diamond blastx -d ${DATABASE} -q ${f} -o matches.tsv -p ${SLURM_CPUS_PER_TASK}
	
#	blastn -num_threads 12 -max_hsps 1 -max_target_seqs 1 -task blastn-short \
#	  -query ${f} \
#	  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen sgi sacc sblastnames staxids sscinames" \
#	  -db ${DATABASE} \
#	  -out ./B_results_blast/${fastafile##*/}_IMGVR_nblast_out.txt -evalue 0.01 -word_size 16	
	
done





