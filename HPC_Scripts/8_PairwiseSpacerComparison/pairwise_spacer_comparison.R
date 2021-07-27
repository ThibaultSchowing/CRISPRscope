#
# Get all spacer cluster 100% id representative pairs and compute the best edit distance between seq-seq, and the complement / reverse complement
# Input: output_file_name, test (bool)
# Author: Thibault Schowing, 2020-2021
#

print("SESSION INFO")
print(sessionInfo())


# Get arguments from 2_results_merger.sh
args <- commandArgs()
#print(args)

input_file_name <- args[6]
output_directory <- args[7]
test_ctrl <- F          # Test for a smaller subset of data

#withr::with_libpaths(new = "./packages", install.packages("tidyverse", repos = "http://cran.us.r-project.org"))
library(tidyverse, lib.loc = "./packages")

library(devtools)
library(stringdist)
library(Biostrings)
library(gtools)
library(tictoc) # quick benchmark
library(data.table) # fwrite


#withr::with_libpaths(new = "./packages", install_github("https://github.com/HenrikBengtsson/future"))
#withr::with_libpaths(new = "./packages", install_github("DavisVaughan/furrr"))


library(future, lib.loc = "./packages")
library(furrr, lib.loc = "./packages")
future::plan(multicore)


library(data.table) # fwrite


print("Packages loaded.")
print(test_ctrl)
print(input_file_name)
print(output_directory)


# Directory cluster side
CRISPRscope_tbl_26 <- readRDS(file = input_file_name)

test_slice <- CRISPRscope_tbl_26 %>% slice_head(n = 5)

#CRISPRscope_tbl_26 <- readRDS(file = "../0_data/pairwise_comparison/SEL26_NCBI_CRISPRscope_tibble.rds")
#CRISPRscope_tbl_26 <- readRDS(file ="../0_data/CRISPRscope_results/SEL26/SEL26_NCBI_CRISPRscope_tibble.rds")

print("Getting unique clusters...")

# Unique cluster ids
unique_cluster_ids <- CRISPRscope_tbl_26 %>% select(cluster_spacer_identity) %>% distinct() %>% unlist()
print("Done!")


# First we build a dataset with the cluster 100% number and the identity sequence (doesn't matter which it is). 
print("Building cluster representative dataset ")
clusters_representatives <- CRISPRscope_tbl_26 %>% 
  group_by(cluster_spacer_identity) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  select(cluster_spacer_identity, SpacerSeq) %>% 
  arrange(cluster_spacer_identity) %>% dplyr::rename(cluster = cluster_spacer_identity, seq = SpacerSeq)

print("Done!")
# We have 16332 clusters. But A - B == B - A so we don't want repetition. expand.grid does the repeat. 
# 16332*(0.5 * 16332) - 16332

print("Getting cluster combinations...")
combinations_clusters <- combn(unique_cluster_ids, 2)

all_pairs_id <- tibble(x = c(combinations_clusters[1,]),
                       y = c(combinations_clusters[2,]))
print("Done!")

print("Check test variable: ")
print(test_ctrl)

if(isTRUE(test_ctrl)){
  print("Test mode: only 500 first rows will be computed. Getting pairs...")
  test_pairs_id <- all_pairs_id %>% slice_head(n = 500)
  
  spacers_representative_pairs <- test_pairs_id %>% 
    left_join(clusters_representatives, by = c("x" = "cluster")) %>% 
    left_join(clusters_representatives, by = c("y" = "cluster"))
  
  # test_pairs %>% write.csv(file = "./OUTPUT/CRISPRscope_tbl_26_pairwise_spacers.csv", col.names = TRUE, sep = ",")
}else{
  print("Normal mode: big computation coming! Getting pairs...")
  spacers_representative_pairs <- all_pairs_id %>%
    left_join(clusters_representatives, by = c("x" = "cluster")) %>% 
    left_join(clusters_representatives, by = c("y" = "cluster"))  
    #write.csv(file = "./OUTPUT/CRISPRscope_tbl_26_pairwise_spacers.csv", col.names = TRUE, sep = ",")
}

print("Done!")



# -------------
# Compute the minimal edit distance for seq - seq | seq - rev(seq) | seq - revComple(seq) - seq 
# -------------


pairwise_comparison <- function(x,y){
  
  # Complement and reverse complement of one of the sequence. 
  
  y.complement <- as.character(complement(DNAStringSet(y)))
  y.reverseComplement <- as.character(reverseComplement(DNAStringSet(y)))
  
  # Returns the minimal distance between the sequence x and the variations of sequence y
  
  return(
    min(stringdist(x,y),
        stringdist(x, y.reverseComplement))
  )
  
}

# Read the csv file containing all the clusters - sequences pairs
# BIG FILE ! (use test file)

# use 1 script option - removed this
# spacers_representative_pairs <- read.csv(file = "./OUTPUT/CRISPRscope_tbl_26_pairwise_spacers.csv", header = TRUE, sep = ",")



#------------------
# Benchmark
#------------------

# Benchmark between map2 and furrr for 150 pairs | 50'000
# Note that performance on windows can be slower on smaller datasets due to thread set up time. 

# tic()
# spacers_representative_pairs %>% mutate(dist = future_map2_dbl(seq.x, seq.y, pairwise_comparison))%>% unnest(dist)
# toc() # 6 sec | 167.7 sec
# 
# tic()
# spacers_representative_pairs%>% mutate(dist = map2_dbl(seq.x, seq.y, pairwise_comparison))%>% unnest(dist)
# toc() # 2.46 sec | 730.4 sec


#----------------
# Quick comparison of pairwise_comparison()
#----------------

# asdf <- spacers_representative_pairs %>% mutate(dist = future_map2_dbl(seq.x, seq.y, pairwise_comparison, .progress = TRUE)) %>% unnest(dist)
# qwer <- spacers_representative_pairs %>% mutate(dist = future_map2_dbl(seq.x, seq.y, pairwise_comparison, .progress = TRUE)) %>% unnest(dist)
# 
# asdf == qwer

#----------------
# Add the distance for all the entries
#----------------
print("Big computation starting! See you later!")
tic()

#fwrite(spacers_representative_pairs %>% mutate(dist = future_map2_dbl(seq.x, seq.y, pairwise_comparison, .progress = TRUE)) %>% unnest(dist), 
#       file="./pairwise_editdistance_spacer_comparison.csv")


splitted_pairs <- split(spacers_representative_pairs, spacers_representative_pairs$y)

conveniant_loop_variable = 1
total_slices = length(splitted_pairs)
for (slice in 1:total_slices) {
  print(paste("[R DEBUG]Slice: ", conveniant_loop_variable, " - Start", sep = ""))
  
  data = splitted_pairs[[slice]]

  

  #print(paste(output_directory,"/slice_",conveniant_loop_variable , "_pairwise_editdistance_spacer_comparison.csv", sep = ""))

  fwrite(data %>%
           mutate(dist = future_map2_dbl(seq.x, seq.y, pairwise_comparison, .progress = TRUE)) %>%
           unnest(dist),
         file=paste("./output_slices/slice_",conveniant_loop_variable , "_pairwise_editdistance_spacer_comparison.csv", sep = ""))

  print(paste("[R DEBUG]Slice: ", conveniant_loop_variable, " - Done", sep = ""))
  conveniant_loop_variable = conveniant_loop_variable + 1
}


print("Big computation done! Congratulations!")
toc()

print("Exit R program.")



