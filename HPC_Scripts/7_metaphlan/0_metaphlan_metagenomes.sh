#!/bin/bash

#SBATCH -o slurm_output/metaphlan-%j-output.txt
#SBATCH -e slurm_output/metaphlan-%j-error.txt

#SBATCH --job-name="mtphlan"
#SBATCH --time=96:00:00
#SBATCH --mem=150G
#SBATCH --partition=pall
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --export=NONE


module add vital-it/7




# singularity run-help /software/singularity/containers/MetaPhlAn2-2.7.8-1.biocontainers.simg;



# singularity exec --bind /data/projects/p539_crisprscope/ /software/singularity/containers/MetaPhlAn2-2.7.8-1.biocontainers.simg metaphlan2.py <options>



BaseLocation_output=/data/projects/p539_crisprscope/9_metaphlan/OUTPUT

#CONTAINER=/software/singularity/containers/MetaPhlAn2-2.7.8-1.biocontainers.simg
#CONTAINER=/data/projects/p539_crisprscope/MetaPhlAn2-2.7.8-1.biocontainers.simg

CONTAINER=/software/singularity/containers/MetaPhlAn-3.0.7-1.ubuntu20.sif
DATE=20210510
rm ${BaseLocation_output}/metaphlan/allStrains_metaphlan.txt

for bioprojectsss in $(ls /data/projects/p539_crisprscope/0_data/metagenomes) #loop through the different bioprojects
do
	echo "======================================="
	echo -e "Bioproject: "${bioprojectsss}
	for biosamplesss in $(ls /data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}) #loop through the different biosamples
	do
		echo "======================================="
		echo -e "Biosample: "${biosamplesss}

		rm -r ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}
		mkdir -p ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}/{bowtie2ouput,metaphlan_profile}

		#check number of read files
		numReadss=$(ls /data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss} |wc -l)
		
		SAMPLE_FOLDER=/data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/		

		if [ "$numReadss" = "2"  ]
		then
			echo -e "there are two read files in the directory, therefore we run a PE metaphlan"
			VAR1="/data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/*_1.fastq.gz"
			VAR2="/data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/*_2.fastq.gz"
			echo ${VAR1}
			echo ${VAR2}
			
			zcat ${VAR1} ${VAR2} > ${SAMPLE_FOLDER}_input.fastq
			INPUT=${SAMPLE_FOLDER}_input.fastq
			
			singularity exec --bind /data/projects/p539_crisprscope/ --bind /db/SOFTWARE/MetaPhlAn/ ${CONTAINER} metaphlan --input_type fastq $INPUT --nproc 5 --bowtie2out  ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}///bowtie2ouput/metagenome_${biosamplesss}_trimmed.bowtie2.bz2 > ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}/metaphlan_profile/profiled_metagenome_${names}_trimmed.txt
			
			rm ${SAMPLE_FOLDER}_input.fastq
			
			grep "s__"  ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}/metaphlan_profile/profiled_metagenome_${names}_trimmed.txt |grep "t__" -v|awk -F "[\t|]" -v culturess="$bioprojectsss" -v strainzzzss="$biosamplesss" '{OFS="\t"}{print culturess,strainzzzss,$7,$15}' >>  ${BaseLocation_output}/metaphlan/${DATE}_allStrains_metaphlan.txt



		else
			echo -e "there are only one read file in the directory, therefore we run a single end bwa alignement"
			
			VAR1="/data/projects/p539_crisprscope/0_data/metagenomes/${bioprojectsss}/${biosamplesss}/*.fastq.gz"
			echo ${VAR1}
			
			zcat ${VAR1} > ${SAMPLE_FOLDER}_input.fastq
			
			INPUT=${SAMPLE_FOLDER}_input.fastq
			
			singularity exec --bind /data/projects/p539_crisprscope/ --bind /db/SOFTWARE/MetaPhlAn/ ${CONTAINER} metaphlan --input_type fastq $INPUT --nproc 5 --bowtie2out  ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}///bowtie2ouput/metagenome_${biosamplesss}_trimmed.bowtie2.bz2 > ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}/metaphlan_profile/profiled_metagenome_${names}_trimmed.txt
			
			rm ${SAMPLE_FOLDER}_input.fastq
			
			grep "s__"  ${BaseLocation_output}/metaphlan/${bioprojectsss}_${biosamplesss}/metaphlan_profile/profiled_metagenome_${names}_trimmed.txt |grep "t__" -v|awk -F "[\t|]" -v culturess="$bioprojectsss" -v strainzzzss="$biosamplesss" '{OFS="\t"}{print culturess,strainzzzss,$7,$15}' >>  ${BaseLocation_output}/metaphlan/${DATE}_allStrains_metaphlan.txt


		fi

	done #biosamplessss
done #bioprojectsss








