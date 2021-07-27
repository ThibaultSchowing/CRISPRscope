#!/bin/bash
#SBATCH --job-name=merge                           # create a short name for your job
#SBATCH --partition=pall
#SBATCH --nodes=1                                    # node count
#SBATCH --ntasks=1                                   # total number of tasks across all nodes
#SBATCH --cpus-per-task=30                           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=20G                            # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00                           # total run time limit (HH:MM:SS)
#SBATCH --error=slurm-output/merge_results.err
#SBATCH --output=slurm-output/merge_results.out

# Merges files from output_slices into one output_file containing 16505 lines of 16505 pairwise comparison edit distance result. 



output_file_30=./merged_slices_custom_30_1.csv
output_file_not30=./merged_slices_custom_not30_2.csv

sorted_output_file_not30=./sorted_results_custom_not30.csv
sorted_output_file_30=./sorted_results_custom_30.csv


rm ${output_file_30}
rm ${output_file_not30}




for slice in $(ls ./output_slices)
do

#	cat ./output_slices/${slice}  | cut -d "," -f 5 | tr '\n' ',' | sed 's/.$/\n/'  | cut -d ',' -f 2-  >> ${output_file}
#	cat ./output_slices/${slice} | tail -n +2 | cut -d "," -f 5 >> ${output_file}

# Instead of just getting the distance, here we also want the length of the two sequences. 

	cat ./output_slices/${slice} | tail -n +2 | awk -F, '{if((length($3) == 30) && length($4) == 30) print }' | cut -d "," -f 5 >> ${output_file_30}
	cat ./output_slices/${slice} | tail -n +2 | awk -F, '{if((length($3) != 30) || length($4) != 30) print }' | cut -d "," -f 5 >> ${output_file_not30}
done


#cat ${output_file} | sort -n | uniq -c > ${sorted_output_file}

cat ${output_file_30} | sort -n | uniq -c > ${sorted_output_file_30}
cat ${output_file_not30} | sort -n | uniq -c > ${sorted_output_file_not30}


