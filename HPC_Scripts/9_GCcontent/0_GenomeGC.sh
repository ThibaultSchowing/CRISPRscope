#!/bin/bash

#SBATCH -o slurm_output/gc-output.txt
#SBATCH -e slurm_output/gc-error.txt

#SBATCH --job-name="GC"
#SBATCH --time=1-00:00:00
#SBATCH --mem=25G
#SBATCH --partition=pall
#SBATCH --export=NONE
#SBATCH --ntasks=1



indir="../0_data/assemblies/MERGED/ANI"

rm -r ./OUTPUT/
mkdir -p ./OUTPUT/

for species in $(ls ${indir})
do

	for genome in $(ls ${indir}/${species}/*.fna)
	do
		s=${genome##*/}
		genome_name=${s%.*}
		#echo $genome
		awk '	BEGIN{RS=">";FS="\n";
			print "name\tA\tC\tG\tT\tN\tlength\tGC%"}
			NR>1{sumA=0;sumT=0;sumC=0;sumG=0;sumN=0;seq="";
			for (i=2;i<=NF;i++) seq=seq""$i; k=length(seq); 
			for (i=1;i<=k;i++) {
				if (substr(seq,i,1)=="T") sumT+=1; 
				else if (substr(seq,i,1)=="A") sumA+=1; 
				else if (substr(seq,i,1)=="G") sumG+=1; 
				else if (substr(seq,i,1)=="C") sumC+=1; 
				else if (substr(seq,i,1)=="N") sumN+=1}; 
			print $1"\t"sumA"\t"sumC"\t"sumG"\t"sumT"\t"sumN"\t"k"\t"(sumC+sumG)/k*100}' ${genome} |tail -n +2 | cut -d '	' -f 1,8 >> ./OUTPUT/${genome_name}_gc_content.tsv
		
		
	done # genome
done # species

