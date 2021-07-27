from Bio import SeqIO
from Bio import pairwise2

import csv
import sys
import collections 
import numpy


print("args: ", sys.argv)
pORG = sys.argv[1]

pMATCH = float(sys.argv[2])
pMISMATCH = float(sys.argv[3])
pOPEN = float(sys.argv[4])
pEXTEND = float(sys.argv[5])
pTHRESHOLD = float(sys.argv[6])

print("Arguments: ", pMATCH, ":", pMISMATCH, ":", pOPEN, ":", pEXTEND, ":", pTHRESHOLD, "-", pORG)





#
#
#
def getCasType(ORG, _MATCH, _MISMATCH, _OPEN, _EXTEND, _PERCENT_THRESHOLD):
	
	
	#_MATCH=1
	#_MISMATCH=-1
	#_OPEN=-.5
	#_EXTEND=0#-.2
	
	#_SCORE_THRESHOLD = 30
	# The score threshold is calculated dynamically as 80% of the shortest sequence's length. 
	#_PERCENT_THRESHOLD = .8
	
	#CCFDBpath="../0_data/DB/CrisprCasFinderDB_dr_1.csv"
	
	# CRISPRmap path (sorry for the too specific name)
	CCFDBpath="../0_data/DB/CRISPRmap/CRISPRmap_data.fa"


	OUTPUT_FILE="./results/" + ORG + "/" + ORG + "_matched_cas_CRISPRmap.csv"
	#ORG="PRJEB6952"
	#ORG="PRJEB15432"
	CRASS_DR="./results/" + ORG + "/" + ORG + ".DR.fasta"
	
	#CCFDB = csv.reader(open(CCFDBpath, "r"), delimiter=";")
	CCFDB=CCFDBpath 
	
	
	output_cas_types = collections.defaultdict(list)
	output_organisms = collections.defaultdict(list)
	output_unmatched = []
	match = False
	
	#
	# Load the dictionary for correspondance
	# Difficult to parse due to "_" being used in different cases (subspecies / groups etc) 
	
	SPECIES_ACRONYM_TAB="../0_data/DB/CRISPRmap/species-acronyms.tab"
	with open(SPECIES_ACRONYM_TAB) as fin:
		rows = ( line.split('\t') for line in fin)
		correspondance_dict = {row[1].strip().split("_")[0]:row[0].strip() for row in rows }
		
			
	#print("Dict access test")
	#for key, value in correspondance_dict.items():
	#	print("-", key, ":", value, ":")
	
	#print("Bdim: ", correspondance_dict.get("Bdim"))
	
	# For each Direct Repeat of the organism:
	for record in SeqIO.parse(CRASS_DR, 'fasta'):
		
#		header=0
		
		#loop through the csv list
		#print("Checking in CCFDB for DR: " + record.seq)
		#for row in CCFDB:
		for row in SeqIO.parse(CCFDB, 'fasta'):
			match = False
			#if the DR matches a DR in the CCFDB
			#No header in fasta
			
#			if(header==1):
#				#skip this row (headers)
#				header=0
#				continue
			
			#print("Row CRISPRmap: ", row.id.split("_")[0], " : ", row.seq)
			#print("CRASS result seq: ", record.id, " : ", record.seq)
			
			# Align the two sequences and get the max score
			# https://biopython.org/docs/1.77/api/Bio.pairwise2.html
			# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec98
			
			#score_only: boolean (default: False). Only get the best score, don’t recover any alignments. 
			#				The return value of the function is the score. Faster and uses less memory.
			
			best_score = pairwise2.align.globalms(record.seq, row.seq, _MATCH, _MISMATCH, _OPEN, _EXTEND, score_only=True)
			
			# Doesn't work as expected, kept for readability. 
			#best_score = pairwise2.align.globalms(record.seq, row[0], match=_MATCH, mismatch=_MISMATCH, open=_OPEN, extend=_EXTEND, score_only=True)
			
			#TODO: check if it makes sense.
			# Threshold: 80% of the smallest/biggest sequence 
			th = int(round(max(len(record.seq), len(row.seq)) * _PERCENT_THRESHOLD ))
			
			#TODO: check for size. Small sequences matching with part of big one will still count as a good match. 
			
			# If the alignment score is higher than the threshold, add (if not already) the matched
			# organism (Cas-type, Superkingdom, Species) to the dictionary (key = DR sequence)
			if(best_score > th):
				 
				 
				# If it is not already present (many times the same cas-type for the same DR)
				if(row.seq + ";" + row.id + ";" + correspondance_dict[row.id.split("_")[0]]  not in output_cas_types[record.seq]):
					output_cas_types[record.seq].append(row.seq + ";" + row.id + ";" + correspondance_dict[row.id.split("_")[0]])
					print(row.seq + ";" + row.id + ";" + correspondance_dict[row.id.split("_")[0]])
					match = True
					 
		
					 
		# Keep track of the unmatched files
		if(not match):
			# record the unmatched sequences somewhere else
			output_unmatched.append(record.seq)
	 
	# Write the DR and their corresponding matche(s) in a csv file.
	print("Writing file!")
	with open(OUTPUT_FILE, 'w') as f:
		f.write("InputDirectRepeat;CRISPRmapDirectRepeat;CRISPRmapID;CompleteID")
		for key in output_cas_types.keys():
			for entry in output_cas_types[key]:
				f.write("\n%s;%s"%(key,entry))
	print("done")














print("start with Bioproject: ", pORG)




getCasType(pORG, pMATCH, pMISMATCH, pOPEN, pEXTEND, pTHRESHOLD)

print("Programm done.. quitting")






