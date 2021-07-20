
library(tidyverse)
library(Biostrings)


#
# Read the organism data (Spacers info)
# Input: species name and source (NCBI or DIALACT Databases)
#
#Lactobacillus_delbrueckii_subsp_bulgaricus_PerSpacer_CRISPR_v1 <- read_spacers_data("Lactobacillus_delbrueckii_subsp_bulgaricus")

args <- commandArgs()
print(args)

spacer_path <- args[6]
repeat_path <- args[7]

print(spacer_path)
print(repeat_path)


#TODO  pass "path" by parameter for each file. 

read_spacers_data <- function(path, source) {
  if (source == "DIALACT") {
    t <-
      read_delim(
        path,
        ";",
        escape_double = FALSE,
        col_names = T,
        trim_ws = TRUE
      ) %>%
      mutate(tmp = ArrayID) %>%
      separate(tmp, c("Strain", "Scaffold", "Array"), sep = "_") %>%
      relocate(Strain, .before = "ShortID") %>%
      relocate(Scaffold, .before = "ShortID") %>%
      relocate(Array, .before = "ShortID") %>%
      relocate(ArrayID, .before = "ShortID") 
  }
  
  if (source == "NCBI") {
    t <-
      read_delim(
        path,  
        ";",
        escape_double = FALSE,
        col_names = T,
        trim_ws = TRUE
      ) %>%
      mutate(tmp = ArrayID) %>%
      separate(tmp, c("Strain", "Scaffold1", "Scaffold2", "Array"), sep = "_") %>%
      unite(Scaffold,
            c("Scaffold1", "Scaffold2"),
            sep = "",
            remove = TRUE,
      ) %>%
      relocate(Strain, .before = "ShortID") %>%
      relocate(Scaffold, .before = "ShortID") %>%
      relocate(Array, .before = "ShortID") %>%
      relocate(ArrayID, .before = "ShortID")
  }
  
  # returns a tibble
return(t)
}

#
# Add to each spacer the relative distance from the leader (or NA if the direction is not set) (0 close to leader, 1 far from leader)
#
add_reldist_leader <- function(input_tibble){
  return(input_tibble %>% 
           add_count(ArrayID, name = "NbSpacerInArray") %>%
           mutate(rel_dist_leader = ifelse(Orientation == "+", SpacerNb/(NbSpacerInArray-1), -1 * (SpacerNb - NbSpacerInArray + 1) / (NbSpacerInArray-1))) %>%
           mutate(rel_dist_leader = ifelse(Orientation == "ND", NA, rel_dist_leader)) %>% 
           mutate(rel_dist_leader = ifelse(NbSpacerInArray == 1, NA, rel_dist_leader)))
}


#
# Converts dna string set to dataframe
#
dss2df <- function(dss) data.frame(length=width(dss), seq=as.character(dss), names=names(dss))


#
# read a fasta file into a tibble
#
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
  
  if (source_type == "NCBI") {
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
# read a fasta file into a tibble
#
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
  
  if (source_type == "NCBI") {
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


















