#!/bin/bash
#SBATCH -o slurm_output/mapping-%j-output.txt
#SBATCH -e slurm_output/mapping-%j-error.txt
#SBATCH --job-name="mapping"
#SBATCH --time=7-00:00:00
#SBATCH --mem=50G
#SBATCH --partition=pall
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --export=NONE
#SBATCH --mail-type=END,FAIL,TIME_LIMIT



module add vital-it/7
module add UHTS/Aligner/bwa/0.7.17;
module add UHTS/Analysis/samtools/1.10;
module add UHTS/Analysis/BEDTools/2.29.2;
module add R/3.6.1;

###-----------------------------------------------------
####Description section
#Here you can define the locations and file names
###-----------------------------------------------------

#mappingFile=<This is the location and file with the DR+spacer+DR>
mappingFile=/data/projects/p539_crisprscope/7_mapping/CRISPRscope_DR_SP_DR.fasta
threads=20

#BaseLocation_output=<directory where your output should go>
BaseLocation_output=/data/projects/p539_crisprscope/7_mapping/OUTPUT

#script_locations=<directory where the R scripts are located>
script_locations=/data/projects/p539_crisprscope/7_mapping/SCRIPTS

###-----------------------------------------------------
####load moduel
#here you load your slurm modules
#the follow tools are needed
###-----------------------------------------------------

##load bwa
##load samtools
##load R ?
##load bedtools

###-----------------------------------------------------
####preperation
#create a bed preperation file
###-----------------------------------------------------

grep ">" ${mappingFile} | sed 's/>//g' |awk -F "\t" '{OFS="\t"}{print $1,"40","60",$1}'  > ${script_locations}/spacer_mapping_file.bed


###-----------------------------------------------------
####script
#here are the indivual scripting steps
###-----------------------------------------------------

##delete and creat output directories
rm -r $BaseLocation_output
mkdir -p $BaseLocation_output/final

#NEW
mkdir -p ${BaseLocation_output}/bwaMapping2DB

##index you mapping file
bwa index $mappingFile



for bioprojectsss in $(ls /data/projects/p539_crisprscope/0_data/metagenomes) #loop through the different bioprojects
do
echo "======================================="
echo -e "Bioproject: "${bioprojectsss}
echo "Bioproject: ${bioprojectsss}" 1>&2
for biosamplesss in $(ls /data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}) #loop through the different biosamples
do
echo "======================================="
echo -e "Biosample: "${biosamplesss}
echo "Biosamplsss ${biosamplesss}" 1>&2

#check number of read files
numReadss=$(ls /data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss} |wc -l)

if [ "$numReadss" = "2"  ]
then
echo -e "there are two read files in the directory, therefore we run a PE bwa alignement"

#TODO: change the *_R2* and *_R1* by *_2.fastq.gz and *_1.fastq.gz 
bwa mem -v 2 -t 37 $mappingFile /data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/*_1.fastq.gz \
/data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/*_2.fastq.gz  | samtools sort -@${threads} -O BAM \
-o ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted.bam - 

else

echo -e "there are only one read file in the directory, therefore we run a single end bwa alignement"
bwa mem -t 37 $mappingFile \
/data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/* | samtools sort -@${threads} -O BAM \
-o ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted.bam - 
fi

##-------------------mapped

samtools view -b -F 4 ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted.bam > \
${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped.bam

##-------------------high quality

samtools view -h -q 30 -b ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped.bam > \
${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30.bam

##-------------------
###-----------------ONLY spacer and not the repeat region mapped reads --> reads mapping to protospacer
##-------------------

rm -r ${BaseLocation_output}/tmp
mkdir -p ${BaseLocation_output}/tmp

cd ${BaseLocation_output}/tmp

##bam header 
samtools view -H ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30.bam  > \
${bioprojectsss}_${biosamplesss}_header.sam  

##creat sam file for calcualtions
samtools view ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30.bam > \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30_sub.sam

## extract reads that map max 40bp
${script_locations}/mappingLengthfiltering_protospacer_thibault.R -i ${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30_sub.sam -o ${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_protospacer.sam

#prepare sam file
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' ${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_protospacer.sam > \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_protospacer_02.sam

#merge header and sam file
cat ${bioprojectsss}_${biosamplesss}_header.sam   \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_protospacer_02.sam |samtools view -S -b > \
${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_protospacer_final.bam


bedtools coverage -bed -a ${script_locations}/spacer_mapping_file.bed -b ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_protospacer_final.bam  >  ${bioprojectsss}_${biosamplesss}_protospacer_coverage.bed

mkdir -p ${BaseLocation_output}/final

awk -F "\t" -v bioprojetzzz="$bioprojectsss" -v biosamplezzz="$biosamplesss" '{OFS="\t"} {print $0,bioprojetzzz,biosamplezzz}' ${bioprojectsss}_${biosamplesss}_protospacer_coverage.bed   > ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_protospacer_coverage.bed

cat ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_protospacer_coverage.bed >> ${BaseLocation_output}/final/allSamples_protospacer_coverage.bed

##-------------------
###-----------------properly mapped reads --> reads mapping to CRISPR SPACER
##-------------------
rm -r ${BaseLocation_output}/tmp
mkdir -p ${BaseLocation_output}/tmp

cd ${BaseLocation_output}/tmp

##bam header 
samtools view -H ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30.bam  > \
${bioprojectsss}_${biosamplesss}_header.sam  

##creat sam file for calcualtions
samtools view ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30.bam > \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30_sub.sam

## extract reads that map max 40bp
${script_locations}/mappingLengthfiltering_spacer_thibault.R -i ${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_mapQ_30_sub.sam -o \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_spacer.sam

#prepare sam file
awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' ${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_spacer.sam > \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_spacer_02.sam

#merge header and sam file
cat ${bioprojectsss}_${biosamplesss}_header.sam   \
${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_spacer_02.sam |samtools view -S -b > \
${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_spacer_final.bam


bedtools coverage -bed -a ${script_locations}/spacer_mapping_file.bed -b \
${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_rawIllumina_bwamem_sorted_mapped_spacer_final.bam  >  \
${bioprojectsss}_${biosamplesss}_spacer_coverage.bed

mkdir -p ${BaseLocation_output}/final

awk -F "\t" -v bioprojetzzz="$bioprojectsss" -v biosamplezzz="$biosamplesss" '{OFS="\t"} {print $0,bioprojetzzz,biosamplezzz}' \
${bioprojectsss}_${biosamplesss}_spacer_coverage.bed   > ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_spacer_coverage.bed

cat ${BaseLocation_output}/bwaMapping2DB/${bioprojectsss}_${biosamplesss}_spacer_coverage.bed >> ${BaseLocation_output}/final/allSamples_spacer_coverage.bed

#-------------------
##-----------------count number of reads for normalization
#-------------------

readss=$(cat /data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/* |grep "^@" -c) 
echo -e ${bioprojectsss}"\t"${biosamplesss}"\t"${readss} >> ${BaseLocation_output}/final/reads4normalization.txt

done #biosamplessss
done #bioprojectsss

