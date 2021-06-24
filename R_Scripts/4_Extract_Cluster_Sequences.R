# Script for Extracting Clustered Gene sequences
library(tidyverse)
library(Biostrings)

e_coli_sequences <- readDNAStringSet("e_coli_peptide.fa")
seq_names <- names(e_coli_sequences)
sequences <- paste(e_coli_sequences)
e_coli_sequences_df <- as_tibble(data.frame(seq_names, sequences))
e_coli_sequences_df[] <- lapply(e_coli_sequences_df, as.character)
e_coli_sequences_df$seq_names <- substr(e_coli_sequences_df$seq_names, 1, nchar(e_coli_sequences_df$seq_names)-2)

fix_target_name <- function(table){
  as_tibble(substr(table$target, 1, nchar(table$target)-2)) # remove last two characters
}

tidy_import <- function(cluster_file){
  cluster_file_tidy <- cluster_file %>%
    select(name) 
  
  cluster_file_tidy$name <- substring(cluster_file_tidy$name, 2)
  cluster_file_tidy$name <- substr(cluster_file_tidy$name, 1, nchar(cluster_file_tidy$name)-1)
  
  return(as_tibble(cluster_file_tidy))
}

join_with_sequences <- function(reference_table, protein_table){
  protein_table %>%
    left_join(reference_table, by = c('name' = 'seq_names'))
}

extract_sequences <- function(clustered_names){
  
  e_coli_sequences <- readAAStringSet("e_coli_peptide.fa") # load from current working directory
  seq_names <- names(e_coli_sequences)
  sequences <- paste(e_coli_sequences)
  e_coli_sequences_df <- as_tibble(data.frame(seq_names, sequences))
  e_coli_sequences_df[] <- lapply(e_coli_sequences_df, as.character)
  e_coli_sequences_df$seq_names <- substr(e_coli_sequences_df$seq_names, 1, nchar(e_coli_sequences_df$seq_names)-2)
  
  HMMER_sequence_hit <- join_with_sequences(e_coli_sequences_df, clustered_names)
  
}

writeFASTAclust <- function(data, filename) {
  fastalines = c()
  for(rowNum in 1:nrow(data)){
    fastalines = c(fastalines, as.character(paste(">", data[rowNum, "name"], sep = "")))
    fastalines = c(fastalines, as.character(data[rowNum, "sequences"]))
  }
  fileConn <- file(filename)
  writeLines(fastalines, fileConn)
  close(fileConn)
}
