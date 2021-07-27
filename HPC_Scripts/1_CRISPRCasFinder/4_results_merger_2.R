print("########################")
print("R merger program 2 start")
print("########################")

library(cli, lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library("seqRFLP", lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(tidyverse , lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(dplyr, lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(Biostrings)

# Get arguments from 2_results_merger.sh
args <- commandArgs()
print(args)

parsed_output_folder <- args[6]
OutputRmerger <- args[7]
filename_output_spacers <- args[8]
filename_output_repeats <- args[9]
DATE <- args[10]
SOURCE_TYPE <- args[11]


print(paste("[DEBUG]: ", sep=""))
print(paste("[DEBUG]: ", sep=""))
print(paste("[DEBUG]: ", sep=""))
print(paste("[DEBUG]: ", sep=""))
print(paste("[DEBUG]: ", sep=""))

folders <- list.dirs(parsed_output_folder, recursive=FALSE)[-1]
print("[DEBUG]: folders")
print(folders)

print("[DEBUG]: init per spacer tibble")

# init per spacer tibble
CRISPRscope_tibble <- tibble()
CRISPRscope_temp_repeat_tibble <- tibble()















print("[DEBUG]: Start looping through organisms")
for (ORGANISM in folders){
	print(paste("[DEBUG]: Organism --> ", ORGANISM[[1]], sep=""))
	# read rds
	orgbasename <- basename(ORGANISM[[1]])
	spacer_tibble <- readRDS(file = paste(ORGANISM,"/", orgbasename, "_spacer_tibble.rds",sep=""))
	repeat_tibble <- readRDS(file = paste(ORGANISM,"/", orgbasename, "_repeat_tibble.rds",sep=""))
	# merge tibbles if not empty
	print(paste("[DEBUG]: tibble dimention spacer", dim(spacer_tibble)[1], sep=""))
	print(paste("[DEBUG]: tibble dimention repeat", dim(repeat_tibble), sep=""))
	
	if(dim(spacer_tibble)[1] > 0){
		CRISPRscope_tibble = bind_rows(CRISPRscope_tibble, spacer_tibble)
		CRISPRscope_temp_repeat_tibble = bind_rows(CRISPRscope_temp_repeat_tibble, repeat_tibble)	
	}
}
print("[DEBUG]: Merging tibbles done")


CRISPRscope_tibble %>% summarise()
print(CRISPRscope_tibble)
print(colnames(CRISPRscope_tibble))
print(CRISPRscope_temp_repeat_tibble)
print(colnames(CRISPRscope_temp_repeat_tibble))

print("[DEBUG]: Joining repeats to main tibble...")
# Merge the repeats to the main tibble 

#
# Join1: join the repeats to the main DF and clean column order
# Next: join the CRISPRCasTyper cas types
#

CRISPRscope_tibble <- left_join(CRISPRscope_tibble, CRISPRscope_temp_repeat_tibble, by = c("Organism", "Strain", "Scaffold", "Array")) %>%
	mutate(DR_length = length, DR_seq = seq, Spacer = SpacerNb, ShortSpacerID = ShortID) %>% 
	select(-length, -seq, -SpacerNb, -ShortID,) %>% 
	relocate(Spacer, .after = Array) %>% 
	relocate(Organism, .before = Strain) %>% 
	relocate(ShortSpacerID, .after = ArrayID) 

print(CRISPRscope_tibble)
print(colnames(CRISPRscope_tibble))
print("[DEBUG]: Joining repeats done.")

print("[DEBUG]: Quality filtering. Removing EvidenceLvl < 4")
# Quality filtering

CRISPRscope_tibble <- CRISPRscope_tibble %>% filter(EvidenceLvl > 3)
print("[DEBUG] Quality filtering done. ")


#
# Export Spacers and Repeats (For CRISPRCasTyper and Clustering)
#
print("[DEBUG]: Start EXPORT of the datasets. ")

# EXPORT FOR CRISPRCasTyper
# 
print("[DEBUG]: Exporting repeats for CRISPRCasTyper")
repeats_unique <- CRISPRscope_tibble %>% select(DR_seq) %>% distinct(DR_seq) 

# TSV used by CRISPRCasTyper
write.table(repeats_unique$DR_seq, paste(OutputRmerger,"/", DATE, "_", SOURCE_TYPE ,"_unique_repeats.tsv", sep=""),na = "", quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
print("[DEBUG]: Done.")

# EXPORT FOR CLUSTERING
print("[DEBUG]: Start fasta EXPORT of the repeats.")
# ------ Repeats
repeats_unique <- repeats_unique %>% rowid_to_column("header")

dfs <- data.frame(repeats_unique$header, repeats_unique$DR_seq)
dfs.fasta = dataframe2fas(dfs, file=paste(OutputRmerger,"/", DATE, "_", SOURCE_TYPE, "_", filename_output_repeats,sep=""))

#TODO: file ???
print("[DEBUG]: Done.")

# ------ Spacers
print("[DEBUG]: Start fasta EXPORT of the spacers.")
spacers_unique <- CRISPRscope_tibble %>% select(SpacerSeq) %>% distinct(SpacerSeq) %>%  
  rowid_to_column("header") %>% add_column(sufx = "genomic") %>% 
  unite(header, c(header, sufx))

dfs <- data.frame(spacers_unique$header, spacers_unique$SpacerSeq)
dfs.fasta = dataframe2fas(dfs, file=paste(OutputRmerger,"/", DATE, "_", SOURCE_TYPE, "_", filename_output_spacers,sep=""))

print("[DEBUG]: Done.")

# Export main dataset
# ------ ALL
print("[DEBUG]: Start EXPORT of the main dataframe.")
saveRDS(CRISPRscope_tibble, file = paste(OutputRmerger,"/", DATE, "_", SOURCE_TYPE, "_CRISPRscope_tibble.rds",sep=""))
print("[DEBUG]: Done.")
print(paste("[DEBUG] : File -> ", paste(OutputRmerger,"/", DATE, "_", SOURCE_TYPE, "_CRISPRscope_tibble.rds",sep=""), sep=""))
























