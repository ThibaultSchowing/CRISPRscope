
#
# Read the organism data (Spacers info)
# Input: species name and source (NCBI or DIALACT Databases)
# Author: Thibault Schowing, 2020-2021
# 

print("SESSION INFO")
print(sessionInfo())
#print("INSTALLED PACKAGES")
#print(installed.packages())

#install.packages("tidyverse", lib="/data/projects/p539_crisprscope/RPACKAGES", repos='http://cran.us.r-project.org')
#install.packages("seqRFLP",  lib="/data/projects/p539_crisprscope/RPACKAGES", repos='http://cran.us.r-project.org')
#install.packages("dplyr",  lib="/data/projects/p539_crisprscope/RPACKAGES", repos='http://cran.us.r-project.org')
#install.packages("cli",  lib="/data/projects/p539_crisprscope/RPACKAGES", repos='http://cran.us.r-project.org')

library(cli, lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library("seqRFLP", lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(tidyverse , lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(dplyr, lib.loc='/data/projects/p539_crisprscope/RPACKAGES')
library(Biostrings)



print("R merger program 1 start")

print("Dplyr version: ")
print(packageVersion('dplyr'))

# Get arguments from 2_results_merger.sh
args <- commandArgs()
#print(args)

spacer_path <- args[6]
repeat_path <- args[7]
source_type <- args[8]
organism_folder <- args[9]
organism <- args[10]

print(spacer_path)
print(repeat_path)
print(source_type)
print(organism_folder)
print(organism)



#-------------------------------------------------------------
#
# Reads the CSV saved by CRISPRCasFinder / parser
# Parameters: - Path of the csv file
# 	      - Type of source (for parsing strings)
#
# Returns a tibble
#-------------------------------------------------------------
read_spacers_data <- function(path, source_type) {
  print("Inside read_spacers_data")
  t <- read_delim(path,";",escape_double = FALSE,col_names = T,trim_ws = TRUE)

  if (source_type == "DIALACT") {
    t <- t %>% mutate(tmp = ArrayID) %>% separate(tmp, c("Strain", "Scaffold", "Array"), sep = "_") %>%
      relocate(Strain, .before = "ShortID") %>%
      relocate(Scaffold, .before = "ShortID") %>%
      relocate(Array, .before = "ShortID") %>%
      relocate(ArrayID, .before = "ShortID") 
  }
  
  if (source_type == "NCBI" || source_type == "TEST") {
	t <- t %>% 
	mutate(tmp = ArrayID) %>% 
	separate(tmp, c("Strain", "Scaffold1", "Scaffold2", "Array"), sep = "_") %>% 
	unite(Scaffold,c("Scaffold1", "Scaffold2"),sep = "",remove = TRUE,) %>% 
	relocate(Strain, .before = "ShortID") %>% 
	relocate(Scaffold, .before = "ShortID") %>% relocate(Array, .before = "ShortID") %>% relocate(ArrayID, .before = "ShortID")
  }
  print("before return")
  return(t)
}

#-------------------------------------------------------------
# Add to each spacer the relative distance from the leader (or NA if the direction is not set) (0 close to leader, 1 far from leader)
#-------------------------------------------------------------
add_reldist_leader <- function(input_tibble){
  print(tbl_sum(input_tibble))
  print(typeof(input_tibble))
  return(input_tibble %>% 
           add_count(ArrayID, name = "NbSpacerInArray") %>%
           mutate(rel_dist_leader = ifelse(Orientation == "+", SpacerNb/(NbSpacerInArray-1), -1 * (SpacerNb - NbSpacerInArray + 1) / (NbSpacerInArray-1))) %>%
           mutate(rel_dist_leader = ifelse(Orientation == "ND", NA, rel_dist_leader)) %>% 
           mutate(rel_dist_leader = ifelse(NbSpacerInArray == 1, NA, rel_dist_leader)))
}

#------------------------------------------------------------- 
# 
# Read a fasta file into a tibble
# Used afterwards to join the repeats to the spacers
#-------------------------------------------------------------
read_DR_fasta_tibble <- function(path, source_type) {
  fastafile = path
  
  dna <- readDNAStringSet(fastafile)
  dfdna <- dss2df(dna)
  tbldna <- as_tibble(dfdna)
  
  print(tbldna)
  
  if (source_type == "DIALACT") {
    print("Source: DIALACT")
    tbldna <- tbldna %>%
      separate(names, c("names", "Array"), sep = '_(?=[0-9]+$)') %>%
      separate(names, c("names", "Scaffold"), sep = '_(?=[[:alnum:]]+$)') %>%
      separate(names, c("names", "Strain"), sep = '_(?=[[:alnum:]-]+$)') %>%
      mutate(Organism = names) %>%
      select(-names) %>%
      relocate(Organism, .before = "Strain") %>%
      relocate(seq, .after = "Array") %>%
      relocate(length, .after = "seq")
  }
  
  if (source_type == "NCBI" || source_type == "TEST") {
    print("Source: NCBI")
    tbldna <- tbldna %>%
      separate(names, c("names", "Array"), sep = '_(?=[0-9]+$)') %>%
      separate(names, c("names", "Scaffold1"), sep = '_(?=[[:alnum:]]+$)') %>%
      separate(names, c("names", "Scaffold2"), sep = '_(?=[[:alnum:]]+$)') %>%
      unite(Scaffold, c("Scaffold2","Scaffold1"), sep = "") %>% 
      separate(names, c("names", "Strain"), sep = '_(?=[[:alnum:].]+$)') %>%
      mutate(Organism = names) %>%
      select(-names) %>%
      relocate(Organism, .before = "Strain") %>%
      relocate(seq, .after = "Array") %>%
      relocate(length, .after = "seq")
    
    
  }
  print("after")
  print(tbldna)
  
  return(tbldna)
}

#
# Converts dna string set to dataframe
#
dss2df <- function(dss) data.frame(length=width(dss), seq=as.character(dss), names=names(dss))



#=============================================================================================#
print("YOOOOO")

print("[DEBUG] --------- Spacers")

spacer_tibble <- add_reldist_leader(read_spacers_data(spacer_path, source_type))

print("[DEBUG] --------- Repeats")
repeat_tibble <- read_DR_fasta_tibble(repeat_path, source_type)



#
# Save the tibble to RDS for each organism
#

print("[DEBUG] : ")
print(paste(organism_folder, "/", organism,"_spacer_tibble.rds", sep=""))
print("[DEBUG]: working dir")
print(getwd())

path_spacer_rds <- paste(organism_folder, "/", organism,"_spacer_tibble.rds", sep="")
saveRDS(spacer_tibble, file = path_spacer_rds)

path_repeat_rds <- paste(organism_folder, "/", organism, "_repeat_tibble.rds", sep="")
saveRDS(repeat_tibble, file = path_repeat_rds)

# ATTENTION: Repeats are exported only after quality filtering and thus only after merging. But we need the repeats from read_DR_fasta_tibble to merge repeats and spacers (they couldn't be parsed together)


# Here every parsed output (in Parsed_output/OrganismX/) has two associated rds files. 
# As this script is executed once for every organism, another script will loop through the outputs and merge the tibbles. 



print("[DEBUG] - END")


























