print("########################")
print("R cluster merger program start")
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
spacer_cluster_csv <- args[7]
repeat_cluster_csv <- args[8]


CRISPRscope_tbl <- readRDS(file=CRISPRscope_main_tibble_path)

Clusters_CRISPRscope_spacers <- read_csv(spacer_cluster_csv) %>% distinct(seq, cluster, identity)
Clusters_CRISPRscope_repeats <- read_csv(repeat_cluster_csv) %>% distinct(seq, cluster, identity)

CRISPRscope_tbl <- CRISPRscope_tbl %>% 
  left_join(Clusters_CRISPRscope_spacers, by = c("SpacerSeq" = "seq")) %>% 
  mutate(cluster_spacer = cluster) %>% 
  mutate(identity_spacer_cluster = identity) %>% 
  select(-cluster, -identity) %>% 
  relocate(cluster_spacer, .after = "SpacerSeq") %>% 
  relocate(identity_spacer_cluster, .after = "cluster_spacer")


CRISPRscope_tbl <- CRISPRscope_tbl %>% left_join(Clusters_CRISPRscope_repeats, by = c("DR_seq" = "seq")) %>% 
  mutate(cluster_repeat = cluster) %>% 
  mutate(identity_repeat_cluster = identity) %>% 
  select(-cluster, -identity) %>% 
  relocate(cluster_repeat, .after = "DR_seq") %>% 
  relocate(identity_repeat_cluster, .after = "cluster_repeat")


print(colnames(CRISPRscope_tbl))

print("Saving rds")
saveRDS(CRISPRscope_tbl, file = CRISPRscope_main_tibble_path)


print("Done joining clusters")
print("exit 7_cluster_merger.R")







