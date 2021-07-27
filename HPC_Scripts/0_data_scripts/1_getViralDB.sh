#
#
#
#
#



date=20201110

#mkdir -p /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt

curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt' > /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt

##------columns-------

#1= assembly_accession	
#2=bioproject	
#3=biosample	
#4=wgs_master	
#5=refseq_category	
#6=taxid	
#7=species_taxid	
#8=organism_name	
#9=infraspecific_name	
#10=isolate	
#11=version_status	
#12=assembly_level	
#13=release_type	
#14=genome_rep	
#15=seq_rel_date	
#16=asm_name	
#17=submitter	
#18=gbrs_paired_asm	
#19=paired_asm_comp	
#20=ftp_path	
#21=excluded_from_refseq	
#22=relation_to_type_material


# Subset


##----------Streptococcus
name="Streptococcus"

grep "${name}" -i /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt | awk -F "\t" 'BEGIN{OFS="\t"}{if($11=="latest" && $12="complete")print $0}' > \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Streptococcus_genomic_file.txt

##----------lactobacillus
name="Lactobacillus"

grep "${name}" -i /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt | awk -F "\t" 'BEGIN{OFS="\t"}{if($11=="latest" && $12="Complete Genome")print $0}' > \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactobacillus_genomic_file.txt

##----------lactococcus-
name="Lactococcus"

grep "${name}" -i /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt | awk -F "\t" 'BEGIN{OFS="\t"}{if($11=="latest" && $12="complete")print $0}' > \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactococcus_genomic_file.txt


##----------Propionibacterium
name="Propionibacterium"

grep "${name}" -i /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt | awk -F "\t" 'BEGIN{OFS="\t"}{if($11=="latest" && $12="complete")print $0}' > \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Propionibacterium_genomic_file.txt

##----------Leuconostoc
name="Leuconostoc"

grep "${name}" -i /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/genomic_file.txt | awk -F "\t" 'BEGIN{OFS="\t"}{if($11=="latest" && $12="complete")print $0}' > \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Leuconostoc_genomic_file.txt



#==========
##-------download fna
##==========

# /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs
# mkdir /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads


#makeblastdb -max_file_sz 10GB -in /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/phage_db.fasta -out /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/phage_db -parse_seqids -dbtype nucl

##-------------download seperate for Ldel 

cat /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactobacillus_genomic_file.txt |cut -f 20  |sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|'> /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactobacillus_download_file.txt

rm -r /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactobacillus
mkdir -p  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactobacillus
cd /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactobacillus


wget --input /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactobacillus_download_file.txt

gunzip /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactobacillus/*fna.gz
cat  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactobacillus/*fna >  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactobacillus_phage_db.fasta


##-------------download seperate for Sterm 

cat /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Streptococcus_genomic_file.txt  |cut -f 20  |sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|'> \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Streptococcus_download_file.txt

rm -r /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Streptococcus
mkdir -p  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Streptococcus
cd /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Streptococcus


wget --input /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Streptococcus_download_file.txt

gunzip /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Streptococcus/*fna.gz
cat  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Streptococcus/*fna >  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Streptococcus_phage_db.fasta


##-------------download seperate for Lactococcus lactis 

cat /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactococcus_genomic_file.txt  |cut -f 20  |sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|'> \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactococcus_download_file.txt

rm -r /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactococcus
mkdir -p  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactococcus
cd /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactococcus


wget --input /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Lactococcus_download_file.txt

gunzip /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactococcus/*fna.gz
cat  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactococcus/*fna >  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Lactococcus_phage_db.fasta



##-------------download seperate for Propionibacterium

cat /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Propionibacterium_genomic_file.txt  |cut -f 20  |sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|'> \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Propionibacterium_download_file.txt

rm -r /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Propionibacterium
mkdir -p  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Propionibacterium
cd /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Propionibacterium


wget --input /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Propionibacterium_download_file.txt

gunzip /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Propionibacterium/*fna.gz
cat  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Propionibacterium/*fna >  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Propionibacterium_phage_db.fasta




##-------------download seperate for Leuconostoc

cat /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Leuconostoc_genomic_file.txt  |cut -f 20  |sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|'> \
/data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Leuconostoc_download_file.txt

rm -r /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Leuconostoc
mkdir -p  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Leuconostoc
cd /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Leuconostoc


wget --input /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/logs/Leuconostoc_download_file.txt

gunzip /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Leuconostoc/*fna.gz
cat  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Leuconostoc/*fna >  /data/projects/p539_crisprscope/0_data/DB/Phage_DB/BLAST/${date}/downloads/Leuconostoc_phage_db.fasta




























