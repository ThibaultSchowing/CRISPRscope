#!/bin/bash

#SBATCH -o makeblastdb-out.txt
#SBATCH -e makeblastdb-err.txt

#SBATCH --job-name="makeDB"
#SBATCH --mem=20G
#SBATCH --partition=pall



module add vital-it/7
module add Blast/ncbi-blast/2.10.1+;



# create a blastable DB from the fasta file

date=20201110
makeblastdb -max_file_sz 3GB -in /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/phage_db.fasta -out /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/phage_db -parse_seqids -dbtype nucl
