#!/bin/bash

#SBATCH --job-name=hmm_search
#SBATCH --partition=pall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/%x_%j.out
#SBATCH --error=slurm_output/%x_%j.err

####--------------------------------------
##preparation
##set you bash variables in order to quickly call them in the script
####--------------------------------------
#username=vsomervi

####---------------------
##modules
####--------------------------------------

#module load HPC/Software 
#module load SequenceAnalysis/HMM-Profile/hmmer/3.1b2

module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2;

module load vital-it/7

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: Test"

echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
#echo -e "Hello " ${username}

HMM_DIR=/data/projects/p539_crisprscope/6_HMMgeneDetection # /users/vsomervi/scratch/vsomervi/SAGE_II/05_Identification_search
FAA_DIR=/data/projects/p539_crisprscope/0_data/assemblies/NCBI/FAA_SEL26/SEL26_PROKKA

rm ${HMM_DIR}/all.scans
rm  ${HMM_DIR}/all.scans_oneline

for genomesss in $(ls ${FAA_DIR}/| grep -v "outgroup")  ##loop for the genomes
do
  
  rm -r ${HMM_DIR}/TMP/${genomesss}/
  mkdir -p ${HMM_DIR}/TMP/${genomesss}/stdouts
  
  echo "=================================================================="
  echo ${genomesss}
  echo "=================================================================="
  
  for hmmsss in $(ls ${HMM_DIR}/HMMs_final/  | grep ".txt$" |sed 's/.txt//g') ##loop for the HMMss
  do
    
    echo "-----------------------------------------------------------------"
    echo ${hmmsss}
    
    #hmmpress /users/vsomervi/scratch/vsomervi/SAGE_II/04_HMMs/${hmmsss}.txt
    hmmsearch --cpu 3 -E 0.001 --tblout ${HMM_DIR}/TMP/${genomesss}/${hmmsss}.scan.tab -o ${HMM_DIR}/TMP/${genomesss}/stdouts/${hmmsss}.scan ${HMM_DIR}/HMMs_final/${hmmsss}.txt ${FAA_DIR}/${genomesss}
    
    grep -v "^#" ${HMM_DIR}/TMP/${genomesss}/${hmmsss}.scan.tab  |awk -v genomzz="$genomesss" -v hmmzz="$hmmsss" '{OFS="\t"}{print genomzz,hmmzz,$3,$6}' >> ${HMM_DIR}/all.scans
    
    #grep -v "^#" /users/vsomervi/scratch/vsomervi/SAGE_II/05_Identification/${genomesss}/${hmmsss}.scan.tab  |awk '{OFS="\t"}{print $3}'| sed -z 's/\n/,/g;s/,$/\n/' |awk -v genomzz="$genomesss" -v hmmzz="$hmmsss" '{OFS="\t"}{print genomzz,hmmzz,$1}'  >> /users/vsomervi/scratch/vsomervi/SAGE_II/05_Identification/all.scans_oneline
    
  done # hmmsss
done #genomessss

