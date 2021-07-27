
#----------------------------------------------------------------------------------------------
# 
# Goal      : Parses the result.json file from CRISPRCasFinder to create exploitable datasets 
#             TODO: maybe clean it a bit as the base file is not used and not complete (faulty start)
#
#
#
# Arguments : Organism name (e.g. Lb_delbruekii)
#
#----------------------------------------------------------------------------------------------



import json
import sys
import os
import re
import traceback

print("Number of arguments:", len(sys.argv), "arguments.")
print("Argument List:", str(sys.argv))


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



organism = sys.argv[1]
indir = sys.argv[2]
outdir = sys.argv[3] + "/"
source = sys.argv[4] # NCBI or DIALACT


print("Organism: " + organism)
print("Output directory: " + outdir)
print("Source: " + source)


# Regex to remove -i1-1 
regex = re.compile(r"-[ip][0-9]-[0-9]")


print("======================================================================")
print("======================================================================\n")

print("Parsing of CCF output for CRISPR Spacers and Cas \n\nOrganism: " + organism)

print("======================================================================")
print("======================================================================")


print("Organism directory: " + indir)





#=====
# CSV
#=====

#header = "Organism; Strain; Array_Id; Evidence_Level; DR_Consensus; Conservation_DRs; Left_FLANK; [spacers]; Right_FLANK"

#output_file_path = outdir + organism + "_CRISPR_v1.csv"

#=====
# Fasta and home made outputs
#=====
output_fasta_path = outdir + organism + "_CRISPR_v1.fasta"
output_short_fasta_path = outdir + organism + "_CRISPR_v1_short.fasta"
output_per_spacer_path = outdir + organism + "_PerSpacer_CRISPR_v1.csv"

# Per genes ???
# here simply Organism + cas type(s)
output_cas = outdir + organism + "_CasType_v1.csv"
output_cas_per_genes = outdir + organism + "_Cas_genes_v1.csv"

output_repeat_fasta_path = outdir + organism + "_DR_v1.fasta"

output_short_repeat_fasta_path = outdir + organism + "_DR_v1_short.fasta"

# Different output strings
fastaString = ""
shortfastaString = ""
spacerString = "Organism;ShortID;EvidenceLvl;Orientation;ArrayID;SpacerNb;SpacerSeq\n"

castypestring = "Organism;Strain;SequenceID;CasType"
casgenesstring = "Organism;Strain;SequenceID;CasType;Gene"

#new
repeatfastastring = ""

shortrepeatfastaString = ""

#print("Removing previous outputs...")
# Clearing before (re)creating the files
#if os.path.exists(output_file_path):
#	os.remove(output_file_path)

print("Done!")

#print("Opening main output file")
#with open(output_file_path, 'w') as outfile:
	#outfile.write(header)
	#print("Headers for main output file: " + header)

print("looping through CCF result directories: " + str(os.listdir(indir)))

for resdir in os.listdir(indir):

	# UNIQUE STRAIN NAME
	# NCBI: GCF.242342354.1
	# DIALACT: FAM23456-i1-1
	if(source == "NCBI"):
		ts = resdir.split(sep='_')
		strain = ts[1] 
	else:
		strain = resdir.split(sep='_')[1]
	
	
	print("\n\n//////////////////////////////////////")
	print("\n\nResult - STRAIN: " + strain)
	print("\n\n//////////////////////////////////////")
	
	print("Try and open file: " + indir + "/" + resdir + "/result.json")
	
	try:
		with open(indir + "/" + resdir + '/result.json', 'r') as f:
			eprint("File f: " + str(f))
			result_dict = json.load(f)
		
			print("File f: " + str(f) + "is loaded")
			
			print("result_dict['Date']: " + result_dict['Date'])
			
			#print(result_dict)
			
			print("Loop through scaffolds: ")
			print("result_dict.get('Sequences'): " + str(result_dict['Sequences']))
			
			for scaffold in result_dict.get('Sequences'):
				print("Scaffold.get(id): " + str(scaffold.get('Id')))
				
				# Something like: FAM24751-i1-1_scf13
				SequenceID = scaffold.get('Id')
				
				print("Getting CAS...")
				for cas in scaffold.get('Cas'):
					print("Cas: " + str(scaffold.get('Cas')))
						
					if(cas != []):
						print("================")
						print(" CAS ")
						print("================")
						
						# First file: just the general cas type 
						castype = cas.get('Type')
						print("Cas geleral type: " + str(castype))
						
						castypestring += "\n" + organism + ";" + strain + ";" + SequenceID + ";" + castype 
						
						# Second file: for all the genes 
						for gene in cas.get('Genes'):
							geneSubType = gene.get('Sub_type')
							print("Gene sub-type: ", geneSubType)
							casgenesstring += "\n" + organism + ";" + strain + ";" + SequenceID + ";" + castype + ";" + geneSubType 
						
						
				print("Getting Crisprs...")		
				print("Scaffold.get('Crisprs'): " + str(scaffold.get('Crisprs')))
						
				for crisprs in scaffold.get('Crisprs'):
					if(crisprs != []):
						print("\n\n")
						print("------------------------------------------------------")
						print("\n")
						
						# UNIQUE ARRAY_ID (strain_scaffold_ID)
						print("Array_id (name): " + crisprs.get('Name'))
						array_id = crisprs.get('Name')
						if(source == "NCBI"):
							array_id = strain + "_" + array_id
						
						# DONE: replace with regular expression so -i1-1 is not the onlyone removed (Illumina sample 1 assembly 1
						# -> For NCBI -> add Strain to IDs
						short_array_id = re.sub(regex,'',array_id)
						if(source == "NCBI"):
							short_array_id = array_id
						
						# EVIDENCE LEVEL
						print("Evidence level: " + str(crisprs.get('Evidence_Level')))
						evidence_lvl = crisprs.get('Evidence_Level')
						
						# CONSERVATION DR
						print(" Conservation DR: " + str(crisprs.get('Conservation_DRs')))
						conservation_dr = crisprs.get('Conservation_DRs')
						
						# CONSENSUS DR
						print("Consensus DR: " + crisprs.get('DR_Consensus'))
						consensus_dr = crisprs.get('DR_Consensus')
						
						# POTENTIAL ORIENTATION
						orientation = crisprs.get('Potential_Orientation')
						print("Orientation: " + str(orientation))
						
						# List of elements of the array
						regions = crisprs.get('Regions')
						
						# TODO: if NCBI or DIALACT 
						if(source == "DIALACT"):
							# Split array_id (Strain-technology-sample-assemblyNb-scaffold)
							print("Split: " + str(array_id.split("_")))
							straintech, scaffold, crispr_array = array_id.split("_")
							strain2, techsample, assemblyNb = straintech.split("-")
						else:
							straintech = ""
							scaffold = ""
							crispr_array = ""
							strain2 = ""
							techsample = ""
							assemblyNb = ""
						
						# Region is a list with all the regions of the crispr array
						print("Loop through regions of the crispr array: ")
						spacers_tab = []
						right_flank = ""
						left_flank = ""
						for region in regions:
							seqType = region.get('Type')
							if(seqType == "Spacer"):
								# SPACER
								print("spacer")
								print(region.get('Sequence'))
								spacers_tab.append(region.get('Sequence'))
							if(seqType == "LeftFLANK"):
								print("Left flank")
								print(region.get('Sequence'))
								left_flank = region.get('Sequence')
							if(seqType == "RightFLANK"):
								print("Right flank")
								print(region.get('Sequence'))
								right_flank = region.get('Sequence')
						
						
						
						
						
						
						# FASTA DR
						
						fasta_dr_entry = ">" + short_array_id + "\n" + consensus_dr + "\n"
						shortrepeatfastaString += fasta_dr_entry
						
						#new
						fasta_dr_entry_long = ">" + organism + "_" + array_id + "\n" + consensus_dr + "\n"
						repeatfastastring += fasta_dr_entry_long
						
						# FASTA
						#
						# For each spacer in "spacers_tab" create a fasta entry. 
						for i in range(len(spacers_tab)):
								
							# TODO: replace with regular expression so -i1-1 is not the onlyone removed (Illumina sample 1 assembly 1
							# FILE_NAME_PATTERN="^FAM[0-9]+[-][ip][0-9][-][0-9].fna$"
							# short_id = array_id.replace('-i1-1','') + "_" + str(i) 
							print("==========REGEX STUFF==========")
							print("Short Array ID: " + short_array_id)
							short_spacer_id = short_array_id + "_" + str(i)
							print("Short spacer ID: " + short_spacer_id)
							print("===============================")
							
							
							
							fasta_entry = ">" + organism + "_" + orientation + "_" + array_id + "_" + str(i) + "\n" + spacers_tab[i] + "\n"
							short_fasta_entry = ">" + short_spacer_id + "\n" + spacers_tab[i] + "\n"
							
							
							
							fastaString += fasta_entry
							shortfastaString += short_fasta_entry
							
							# Other funny file / fasta header without > and with the sequence on the same line
							# IMPORTANT: the short entry can serv as unique identifier and will be used to produce the v2 of the file -> with clusters
							# IMPORTANT: don't forget to update header
							spacerEntry =  organism + ";" + short_spacer_id + ";" + str(evidence_lvl) + ";" + orientation + ";" + array_id + ";" + str(i) + ";" + spacers_tab[i] + "\n" 
							spacerString += spacerEntry
	except Exception as e:
		print(e)
		
		traceback.print_exc()
		print("CRISPRCasFinder json result file returning an error for strain " + strain + " ... ignoring.")
		print("Switching to next file.")
		continue



print("--------------------------")
print("Parsing done!")
print("Writing files...")

with open(output_fasta_path, 'w') as fastaOutput:
	fastaOutput.write(fastaString)
							
with open(output_short_fasta_path, 'w') as shortfastaOutput:
	shortfastaOutput.write(shortfastaString)

with open(output_per_spacer_path, 'w') as perSpacerOutput:
	perSpacerOutput.write(spacerString)


with open(output_short_repeat_fasta_path, 'w') as shortDRfastaOutput:
	shortDRfastaOutput.write(shortrepeatfastaString)

#new
with open(output_repeat_fasta_path, 'w') as DRfastaOutput:
	DRfastaOutput.write(repeatfastastring)

# CAS
with open(output_cas, 'w') as casoutput:
	casoutput.write(castypestring)

with open(output_cas_per_genes, 'w') as casgenoutput:
	casgenoutput.write(casgenesstring)

print("Done! End of parser script.")


