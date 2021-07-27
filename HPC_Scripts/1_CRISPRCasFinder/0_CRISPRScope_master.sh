#!/bin/bash

#SBATCH -o slurm_output/master-ccf-%j-output.txt
#SBATCH -e slurm_output/master-ccf-%j-error.txt

#SBATCH --job-name="Xtract"
#SBATCH --time=4-00:00:00
#SBATCH --mem=25G
#SBATCH --partition=pall
#SBATCH --export=NONE
#SBATCH --ntasks=1



# Load modules

module load vital-it/7




# Please check https://hpc.nih.gov/docs/job_dependencies.html for pipelining with sbatch

#TODO: verify parameters

ORGANISM=$1



#       activate the environment before executing python scripts
# -> obviously export the .yml or else file before to create a reproducible env. 



#=====================================================
#   Create conda environment #TODO move that in INIT.sh when done
#=====================================================

# check: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments

# anaconda must be installed (conda) 

#source activate pycrispr

#conda deactivate




#=====================================================
# INIT
#=====================================================
#below

#=====================================================
# CCF
#=====================================================

echo "========================================================================="
echo "Start crisprcasfinder"
echo "========================================================================="
#sleep 2


#CCF1=$(sbatch --wait 1.2_CRISPRCasFinder_task_launcher.sh ${ORGANISM})
#echo ${CCF1}
#CCF1_JID=${CCF1##* }



if [ -z "$1" ]
then
        echo "No argument supplied - Give the name of an organism existing in ${SOURCE}"
        exit
fi


ORGANISM=$1


source CRISPRscope.conf


#
# SKIP CCF IF just parsing needed
#
#parsonly=false
if test "$parsonly" = false; then


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Sequence names formating (removes spaces)
#
# Keep only first word for sequences names. CRISPRCasFinder Needs NO SPACES in the sequences names.
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------



echo "---------------------------"
echo "Sequence names verification"
echo "---------------------------"

sed 's/\s.*$//' -i ${INDIR}/*.fna

echo "... Done ! "
echo " "


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Output managment (USE SNAKEMAKE AND REMOVE !)
#
# Removes result if "force_new" is set to true and recompute all
# Creates Output directory if needed
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


if [[ -d "${OUTDIR}" && "$force_new" = true  ]] ; then
        echo "-----"
        echo "Output directory already exists..."
        echo "Deleting all results...."
        rm -rf ${OUTDIR}
        echo "Done!"
        echo "Creates output directory ${OUTDIR}"
        mkdir -p ${OUTDIR}
        echo "Done!"
        echo "-----"
elif [[ -d "${OUTDIR}" && "$force_new" = false ]] ; then
        echo "-----"
        echo "Output directory already exists..."
        echo "Keeping old results..."
        echo "No need to create directory, verify later for reusults..."
        echo "-----"
elif [[ ! -d "${OUTDIR}" ]] ; then
        echo "-----"
        echo "Creating output directory..."
        mkdir -p ${OUTDIR}
        echo "Done!"
        echo "-----"
fi


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# Create EXECUTION directory for TSKS script. Will create a sbatch for the requested organism
# and will execute it in its own folder
#
#
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

# Execution directory

EXEC_DIR=EXECDIR_${SOURCE_TYPE}_${ORGANISM}
TSKS_FILE=tsks_${SOURCE_TYPE}_${ORGANISM}.sh

rm -r ${EXEC_DIR}
mkdir -p ${EXEC_DIR}


echo "#!/bin/bash" > ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH -o ../slurm_output/tsks_${ORGANISM}-%j-output.txt" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH -e ../slurm_output/tsks_${ORGANISM}-%j-error.txt" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --job-name=\"TSKS\"" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --time=72:00:00" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --mem=25G" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --partition=pall" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --cpus-per-task=4" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --nodes=1" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "#SBATCH --export=NONE" >> ${EXEC_DIR}/${TSKS_FILE}
#echo "module add vital-it/7" >> ${EXEC_DIR}/${TSKS_FILE}


for f in $(ls ${INDIR}/*.fna)
do

        filename="${f##*/}"
        filenameWithoutExtension="${filename%.*}"
        echo "Preparing for file ${f}"
        echo "Filename without extension: $filenameWithoutExtension"

        #---------------------------------------
        # Check if results are already present
        #---------------------------------------

        result_missing=true

        for resdir in $(ls ${OUTDIR})
        do
                # Check if the strain result is already there: If the force_new had already been set to true, the directory should be empty.
                #                                              If the strain has already been analysed, it will ignore the file.
                if [[ $resdir == *"$filenameWithoutExtension"* ]];then
                        echo "----------------------"
                        echo "Result already present: "
                        echo "$filenameWithoutExtension was already analysed. Results in $resdir"
                        echo "----------------------"
                        result_missing=false
                        break
                fi
        done

        # IMPORTANT: The -out OUTDIR must be unique for each assembly / fasta file

        SEQOUTDIR=../Output_CRISPRCasFinder/${DATE}/${SOURCE_TYPE}/${ORGANISM}/Result_${filenameWithoutExtension}

        echo "Sequence specific output directory: ../${SEQOUTDIR}"

        if [[ "$result_missing" = true ]]; then
                echo "Results for $filenameWithoutExtension need to be calculated... "
                echo "Adding file to SBATCH file."
                echo "singularity exec -B /data/projects/p539_crisprscope  ${CONTAINER} CRISPRCasFinder.pl -so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -out ${SEQOUTDIR} -in ../$f" >> ${EXEC_DIR}/${TSKS_FILE}

        fi

done


if test "$launch" = true; then
        echo "==============================="
        echo " Launching CRISPRCasFinder for $ORGANISM"
        echo "==============================="
        cd ${EXEC_DIR}
        # The wait should prevent the task_launcher from returning (!dependency!)
        #CCF=$(sbatch --wait ${TSKS_FILE})
	#CCF_JID=${CCF##* }
	#echo ${CCF}
        #sleep 1
	chmod 755 ${TSKS_FILE}
	./${TSKS_FILE}
	
        cd ..

fi


#--------
echo " "
echo "crisprcasfinder done"
echo "-------------------------------------------------------------------------"
echo ". "
echo ". "





fi # exec parsing only 

#=====================================================
# PARSING
#=====================================================


# The python parser uses the FAMXXXXX_i1-1.fasta name format to extract information. It would be a good idea to check for this format and if the format is not matching to add a dummy x0-0 in the end of each file/sequence. 
# Instead of reading the variables (paths, organism, outdir, etc) from the config file, just pass them by CLI



# IMPORTANT: might be necessary to execute the pythons cript from a SBATCh to get the dependency


echo "========================================================================="
echo "Start parsing"
echo "========================================================================="
sleep 3

# problem with --dependency=afterany:${CCF1_JID}
#echo "Launches parsing job: dependency on ${CCF}"
#PARS1=$(sbatch --wait --dependency=afterany:${CCF} 2_0_parser.sh ${ORGANISM})  
#PARS1_JID=${PARS1##* }
#echo ${PARS1_JID}

source activate pycrispr
python3 ./1_CCF_Parser.py ${ORGANISM} ${OUTDIR} ${PARSDIR} ${SOURCE_TYPE}
conda deactivate


echo "parsing done"
echo "-------------------------------------------------------------------------"
echo ". "
echo ". "

#=====================================================
# CLEANING
#=====================================================

# NOT NOW -> Remove CRISPRCasFinder outputs after parsing ? (would need to be computed again..)


# Remove EXECDIR if the program was executed
if test "$launch" = true; then
	rm -r ${EXEC_DIR}
fi








echo "CRISPRCasFinder and Parsing done. Check output in ParsedOutput directory. "

# As this script is executed for ONE organism, the merging has to be executed for all the others
# especially as the CLUSTERING has to be done on all spacers and repeats together. 

# NEXT STEPS 

#=====================================================
# R script V1 -> For each organism Get main dataset and Repeats
#=====================================================

# MERGE HERE to have one Repeat dataset only

#=====================================================
# Using conda get Cas subtype with Repeattyper
#=====================================================


#=====================================================
# 
#=====================================================


# CLUSTERING HAS TO BE DONE WITH ALL THE SEQUENCES. IT IS THEN DONE SEPARATELY WHEN NEEDED









#EOF
