#!/bin/bash
#SBATCH --job-name=paircmp                           # create a short name for your job
#SBATCH --partition=pall
#SBATCH --nodes=1                                    # node count
#SBATCH --ntasks=1                                   # total number of tasks across all nodes
#SBATCH --cpus-per-task=38                           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=80G                            # memory per cpu-core (4G is default)
#SBATCH --time=28-00:00:00                           # total run time limit (HH:MM:SS)
#SBATCH --error=slurm-output/pairwise_comparison.err
#SBATCH --output=slurm-output/pairwise_comparison.out


# IMPORTANT
# When not testing set partition/cpus/mem/time accordingly

module add vital-it/7
module add R/3.6.1;

# Copy data locally on SCRATCH partition

input_file=/data/projects/p539_crisprscope/0_data/pairwise_comparison/export_complete_SEL26.rds
local_input_file=${SCRATCH}/export_complete_SEL26.rds

output_directory="./output_slices"
mkdir -p ${output_directory}

echo "[DEBUG] Copying data to scratch"
# Copy the tibble to the SCRATCH partition
cp ${input_file} ${local_input_file} 


start=`date +%s`
echo "[DEBUG] Launching R script"



Rscript pairwise_spacer_comparison.R ${local_input_file} ${output_directory}




echo "[DEBUG] R script done"
end=`date +%s`
runtime=$((end-start))

echo "[DEBUG] Time: " ${runtime}
