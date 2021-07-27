print("########################")
print("R CCTYPER merger program start")
print("########################")

library(cli, lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library("seqRFLP", lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(tidyverse , lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(dplyr, lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
#library(Biostrings)

# Get arguments from 2_results_merger.sh
args <- commandArgs()
print("--------------")
print(args)
print("--------------")

CRISPRscope_main_tibble_path <- args[6]
cctyper_output_tsv <- args[7]

print("[DEBUG] Loading rds")
CRISPRscope_tbl <- readRDS(file=CRISPRscope_main_tibble_path)




CRISPRscope_cctyper_cas_types <- read_table2(cctyper_output_tsv, col_names=c("Sequence","Cas_subtype","Cas_prediction_probability")) 	

print(CRISPRscope_cctyper_cas_types)
	
#  mutate(Cas_prediction_probability = Prediction, Cas_subtype = Subtype) %>% 
#  select(-probability, -Prediction, -Subtype)


#
# Merging
#
CRISPRscope_tbl <- CRISPRscope_tbl %>% left_join(CRISPRscope_cctyper_cas_types, by = c("DR_seq" = "Sequence"))
print("[DEBUG] - main dataframe after merging with inferred cas-subtypes")
print(CRISPRscope_tbl)



print("[DEBUG] Saving rds")
saveRDS(CRISPRscope_tbl, file = CRISPRscope_main_tibble_path)

print("END")
print("======================")








