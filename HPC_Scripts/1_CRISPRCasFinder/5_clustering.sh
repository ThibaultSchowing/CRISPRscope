#! /bin/bash


#SBATCH -o slurm_output/cluster-%j-output.txt
#SBATCH -e slurm_output/cluster-%j-error.txt

#SBATCH --job-name="clstr"
#SBATCH --time=00:10:00
#SBATCH --mem=25G
#SBATCH --partition=pshort
#SBATCH --export=NONE
#SBATCH --ntasks=1

module add vital-it/7
module add UHTS/Analysis/cd-hit/4.6.8;

input_fasta=$1
OutputRmerger=$2

# TODO: verify input (file type, number etc)

#filename=${input_fasta%.*}
s=${input_fasta##*/}
filename=${s%.*}
extension="${input_fasta##*.}"

INPUT=${input_fasta}

# Create output directory
OUTDIR_BASE=./OUTPUT_CLUSTER_TMP # to delete in the end
OUTDIR=./OUTPUT_CLUSTER_TMP/${filename}


echo "------------------------------"
echo "CLUSTERING with CD-HIT-EST"
echo "------------------------------"

echo "."
echo "."
echo "Input argument (must be a fasta file): ${input_fasta}"
echo "Filename: ${filename}"
echo "Extension: ${extension}"
echo "Input: ${INPUT}"
echo "Output directory: ${OUTDIR}"
echo "."
echo "------------------------------"


echo "START"
echo "Remove / create output directory..."

rm -r ${OUTDIR}
mkdir -p ${OUTDIR}

echo "Done!"


# Create directory for singleton clusters (easier to parse)
echo "Create director for singleton clusters..."
mkdir -p ${OUTDIR}/singleton_clusters
echo "Done!"


# CD-HIT-EST

echo "Start cd-hit-est..."
cd-hit-est -d 0 -i ${INPUT} -o ${OUTDIR}/CLUSTER_${filename}.out -c 0.8
echo "Done!"



# For parsing, separate the cluster file into multiple files
echo "Split cd-hit-est output into multiple files..."
CLST=${OUTDIR}/CLUSTER_${filename}.out.clstr

awk -v a="$OUTDIR" '/>Cluster/{x=""'"a"'"/singleton_clusters/F"++i-1;}{print > x;}' $CLST

echo "Done!"
echo "."
echo "."
echo "Clustering done ! "
echo "-----------------------------"

echo "###########"
echo "Launching python parser..."

# parse output to generate sequence-cluster dictionary



source activate pycrispr
python3 ./6_cluster_parser.py ${input_fasta} ${filename} ${OutputRmerger}
conda deactivate

echo "."
echo "."
echo "."
echo "Parsing done !"
echo "Output file in ${OutputRmerger}"
echo "###########"


#TODO Remove $OUTDIR as the output file is in OutputRmerger 

echo "Remove output directory ${OUTDIR}"
rm -r ${OUTDIR_BASE}
echo "Done!"

echo "Exit clustering.sh"

















