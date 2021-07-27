#
# Converts cd-hit-est results to csv. Import in R in CRISPRscope
#
# 
#
# Arguments: filename -> fasta file. 
#



import pandas as pd
import os
import sys
from Bio import SeqIO

print("=====================")
print("Python cluster parser")
print("=====================")
print("")

print('#Args: ' + str(len(sys.argv)))
print('Arg list: ' + str(sys.argv))


filename=sys.argv[1]
print("Filename: " + filename)

filebasename = sys.argv[2]
print("File base name: " + filebasename)

OutputRmerger = sys.argv[3]
print("output R merger: " + OutputRmerger)

INDIR="./INPUT"
OUTDIR="./OUTPUT_CLUSTER_TMP"


#
# Step 1: read the fasta file to have a Header - Sequence dict
#

with open(filename) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seqlist = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
                identifiers.append(seq_record.id)
                seqlist.append(str(seq_record.seq))

# converts the two lists to a dataframe

d = {'seqid': identifiers, 'seq': seqlist}
df = pd.DataFrame(d)

# Dataframe containing ID - SEQ
print(df.head())


#
# Loop through all the cluster files (Fx) and get Seq. and Cluster IDs
#

directory = OUTDIR + "/" + filebasename + "/singleton_clusters/" 

cdhit_seqID = []
cdhit_clusterID = []
cdhit_identity = []


# For all the clusters (files named FXXX)
for fn in os.listdir(directory):
	# for each clustering file, parse the file and get the DR id and the cluster name

	# For every line
		# if Cluster x -> x in df
	# Separate line
	# print("Working on " + fn)
	with open(directory + fn, "r") as f:
		line = f.readline()
		cnt=1
		while line:
			# do things here
			if(cnt==1):
				# cluster name line
				clusterID=line.split(" ")[1].replace('\n','').strip()
			else:
				# is_ref: 1 if reference genome, similarity % otherwise
				similarity=line.split(" ")[2].replace('\n','').strip()
				is_ref=0

				if(similarity=="*"):
					is_ref=1
				else:
					# When it is not * the split counts the 'at' as a column
					similarity=line.split(" ")[3].replace('\n','').strip()
					# If it is not the reference genome, possible to parse the percentage of identity
					# at -/84.78%
					s = float(similarity[2:-1])
					if(similarity[0] == '-'):
						s*=-1

					is_ref=s

				seqID=line.split(" ")[1][1:-3]

				#df.loc[df['id'] == seqID, 'cluster'] = clusterID
				#df.loc[df['id'] == seqID, 'is_ref'] = is_ref

				#usvdf.loc[usvdf['ShortID'] == seqID, 'cluster'] = clusterID
				#usvdf.loc[usvdf['ShortID'] == seqID, 'is_ref'] = is_ref
				
				#
				# Append entries to the lists
				#
				#print("Sequence ID:" + seqID)
				cdhit_seqID.append(seqID)
				#print("Cluster ID: " + clusterID)
				cdhit_clusterID.append(clusterID)
				#print("Reference: " + str(is_ref))
				cdhit_identity.append(is_ref)
			
			line = f.readline()
			cnt+=1
#
# Create dict and df
#
cdhit_dict = {'seqid': cdhit_seqID, 'cluster': cdhit_clusterID, 'identity':cdhit_identity}
cdhit_df = pd.DataFrame(cdhit_dict)

# Merge with sequences
# df        -> seqid  - sequence
# cdhit_df  -> seqid  - fluster  - identity

print("Fasta dataframe (seqid - seq)")
print(df.head())
print("Cluster dataframe (seqid, clusterid, ref)")
print(cdhit_df.head())

outputdf = pd.merge(df, cdhit_df, how='left', on='seqid')

print("==================")
print("Output dataframe: ")
print(outputdf.head())

print("Writing file...")

outputdf.to_csv(OutputRmerger + "/Clusters_" + filebasename + ".csv", sep = ',', header = True, index = False)


print("...")
print("DONE!")







































