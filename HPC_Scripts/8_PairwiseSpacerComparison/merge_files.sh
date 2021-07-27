#!/bin/bash
#SBATCH --job-name=merge                           # create a short name for your job
#SBATCH --partition=pall
#SBATCH --nodes=1                                    # node count
#SBATCH --ntasks=1                                   # total number of tasks across all nodes
#SBATCH --cpus-per-task=38                           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=20G                            # memory per cpu-core (4G is default)
#SBATCH --time=12:00:00                           # total run time limit (HH:MM:SS)
#SBATCH --error=slurm-output/merge_results.err
#SBATCH --output=slurm-output/merge_results.out

# Merges files from output_slices into one output_file containing 16505 lines of 16505 pairwise comparison edit distance result. 



output_file=./merged_slices.csv
sorted_output_file=./sorted_results.csv


rm ${output_file}


for slice in $(ls ./output_slices)
do

#	cat ./output_slices/${slice}  | cut -d "," -f 5 | tr '\n' ',' | sed 's/.$/\n/'  | cut -d ',' -f 2-  >> ${output_file}
	cat ./output_slices/${slice} | tail -n +2 | cut -d "," -f 5 >> ${output_file}

done


cat ${output_file} | sort -n | uniq -c > ${sorted_output_file}



