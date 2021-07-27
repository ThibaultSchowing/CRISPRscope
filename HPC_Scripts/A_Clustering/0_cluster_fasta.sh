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


# TODO: verify input (file type, number etc)

filename=${input_fasta%.*}
extension="${input_fasta##*.}"

INPUT=./INPUT/$input_fasta

# Create output directory

OUTDIR=./OUTPUT/${filename}

IDENTITY=1

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
cd-hit-est -d 0 -s 0.9 -n 5 -i ${INPUT} -o ${OUTDIR}/CLUSTER_${filename}.out -c ${IDENTITY}
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

#TODO: parse output to generate sequence-cluster dictionary



source activate pycrispr
python3 1_cd-hit-est_parser.py ${input_fasta} ${filename} ${IDENTITY}
conda deactivate

echo "."
echo "."
echo "."
echo "Parsing done !"
echo "Bye!"
echo "###########"





















