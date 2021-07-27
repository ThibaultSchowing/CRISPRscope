#!/bin/bash

#SBATCH -o slurm_output/fastani-output.txt
#SBATCH -e slurm_output/fastani-error.txt

#SBATCH --job-name="fastANI"
#SBATCH --time=7-00:00:00
#SBATCH --partition=pall
#SBATCH --export=NONE
#SBATCH --nodes=1           ### Number of Nodes
#SBATCH --mem=0
#SBATCH --ntasks-per-node=1 ### Number of tasks (MPI processes)
#SBATCH --cpus-per-task=8  ### Number of threads per task (OMP threads


# Load modules

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

module add vital-it/7
module add UHTS/Analysis/fastANI/1.1;

mkdir -p slurm_output
##fastaANI


###---------------------------
##make a list of locations of all Sterm and ldel genomes to comapre
###---------------------------

##-------------------------------all


rm -r /data/projects/p539_crisprscope/8_fastani/{MERGED_run,log}

mkdir -p /data/projects/p539_crisprscope/8_fastani/{MERGED_run,log}

DBss=MERGED
datesss=ANI
for speciesss in $(ls /data/projects/p539_crisprscope/0_data/assemblies/${DBss}/${datesss})
do
	rm /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani_${speciesss}.txt
	echo "===================================="
	echo ${speciesss}
   	ls /data/projects/p539_crisprscope/0_data/assemblies/${DBss}/${datesss}/${speciesss} |wc -l
	ls /data/projects/p539_crisprscope/0_data/assemblies/${DBss}/${datesss}/${speciesss} |grep -c ".fna"
    	for genomessss in $(ls /data/projects/p539_crisprscope/0_data/assemblies/${DBss}/${datesss}/${speciesss})
    	do
      
		genome_short=$(echo ${genomessss}|sed 's/.fna//g')
		echo -en "/data/projects/p539_crisprscope/0_data/assemblies/"${DBss}"/"${datesss}"/"${speciesss}"/"${genome_short}".fna\n"  >> /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani_${speciesss}.txt
		echo -en ${speciesss}"\t"${genome_short}  >> /data/projects/p539_crisprscope/8_fastani/log/filenames_all_description_genomes.txt
      
	done
	
	wc -l /data/projects/p539_crisprscope/8_fastani/log/filenames_all_description_genomes.txt
	wc -l /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani_${speciesss}.txt
	
	# launch here to do per-species fastANI
	srun --export=ALL fastANI --fragLen 1000 -t 8 --ql /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani_${speciesss}.txt \
	--rl /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani_${speciesss}.txt \
	-o /data/projects/p539_crisprscope/8_fastani/MERGED_run/MERGED_ANI_output_${speciesss}.txt
	
	
done

echo "========================"

#wc -l /data/projects/p539_crisprscope/8_fastani/log/filenames_all_description_genomes.txt
#wc -l /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani.txt


###---------------------------
##fastaANI
###---------------------------
#run --export=ALL fastANI --fragLen 1000 -t 8 --ql /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani.txt \
# --rl /data/projects/p539_crisprscope/8_fastani/log/filenames_all_fastani.txt \
# -o /data/projects/p539_crisprscope/8_fastani/20210511_run/20210511_ANI_output_all.txt




