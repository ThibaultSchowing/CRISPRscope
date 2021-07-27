#!/bin/bash

#SBATCH -o slurm_output_nt/blast-output.txt
#SBATCH -e slurm_output_nt/blast-error.txt
#SBATCH --time=12-00:00:00
#SBATCH --job-name="BLASTnt"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=37
#SBATCH --mem=20G
#SBATCH --partition=pall


# Blast the spacers against the NCBI nt database 



module add vital-it/7

module rm Blast/ncbi-blast/2.10.1+;
module add Blast/ncbi-blast/2.9.0+;

module rm Blast/ncbi-blast/2.7.1+;


# makeblastdb -max_file_sz 1GB -in xxxxxxx.fasta -out yyyyyyyy.fasta




#Directory where to copy fasta files of interest
DATA_DIR="./A_splited_cluster_spacers"

#DATABASE=../0_data/DB/PhageWeb/PhageWebDB
#DATABASE=../0_data/DB/IMGVR/IMG_VR_2020-10-12_5.1/IMGVR
#DATABASE=/data/databases/ncbi-blastdbs/nt_2021_04_01
DATABASE=../0_data/DB/NT/nt

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

mkdir -p ./E_results_blast_nt
rm ./E_results_blast_nt/*
#touch ./E_result_blast_nt/blastn_results.txt

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
          -out ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt -evalue 10 -word_size 7
	
	cat ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt >> ./TMP_NCBIntcollection_concat_blastn_results.txt 
	rm ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt
	
	# if empty
	#if [ ! -s ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt ] ; then
	#	rm ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt
	#fi
done

#echo "parallel ?"
#
#ls ${DATADIR} | parallel -a - blastn -num_threads 35 -max_hsps 1 -max_target_seqs 1 -task megablast -show_gis \
#  -query ${mappingFile} \
#  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen sgi sacc sblastnames staxids sscinames" \
#  -db ${DATABASE} \
#  --out ./E_results_blast_nt/{.}.out \
#  -evalue 0.01 -word_size 16
#
#ls ${DATADIR} | parallel -a - blastn -query {} -db ${DATABASE} -num_threads 5 -max_hsps 1 -max_target_seqs 1 -task blastn-short -show_gis -out ./E_results_blast_nt/{.}_NCBIntcollection_nblast_out.txt -evalue 10 -word_size 7
#"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen sgi sacc sblastnames staxids sscinames" -out ./E_results_blast_nt/{.}_NCBIntcollection_nblast_out.txt -evalue 10 -word_size 7

echo "maybe..."

#cat ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt >> ./TMP_NCBIntcollection_concat_blastn_results.txt
#rm ./E_results_blast_nt/${fastafile##*/}_NCBIntcollection_nblast_out.txt


# Merge results later







# Include sscinames etc in the request
# -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle' -query 

#TODO: check previous blast for vincent's output


