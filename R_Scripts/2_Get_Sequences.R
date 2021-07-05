library(tidyverse)
library(Biostrings)
library(ggdendro)
library(reshape)
library(grid)
library(ggthemes)

# WHAT DOES THIS SCRIPT DO? 
# -- CONVERTS THE PROTEIN TABLES TO FASTA SEQUENCES TO ALLOW PARSING -- 

source("scripts/functions_script_ecoli_efflux_diversity_analysis.R")

# extract sequences from HMMER hits 
ABC_sequence_hits <- extract_sequences_from_e_coli(ABC_proteins_full_clean2)
ACR_sequence_hits <- extract_sequences_from_e_coli(ACR_proteins_cut_clean2)
MatE_sequence_hits <- extract_sequences_from_e_coli(MatE_proteins_cut_clean2)
MFS_sequence_hits <- extract_sequences_from_e_coli(MFS_proteins_cut_clean2)
SMR_sequence_hits <- extract_sequences_from_e_coli(SMR_proteins_full_clean2)

# convert the second column to a character 
ABC_sequence_hits[, 2] <- apply(ABC_sequence_hits[, 2], 2, function(x) as.character(as.factor(x)))
ACR_sequence_hits[, 2] <- apply(ACR_sequence_hits[, 2], 2, function(x) as.character(as.factor(x)))
MatE_sequence_hits[, 2] <- apply(MatE_sequence_hits[, 2], 2, function(x) as.character(as.factor(x)))
MFS_sequence_hits[, 2] <- apply(MFS_sequence_hits[, 2], 2, function(x) as.character(as.factor(x)))
SMR_sequence_hits[, 2] <- apply(SMR_sequence_hits[, 2], 2, function(x) as.character(as.factor(x)))

# Convert df to FASTA
writeFASTA <- function(data, filename) {
  fastalines = c()
  for(rowNum in 1:nrow(data)){
    fastalines = c(fastalines, as.character(paste(">", data[rowNum, "target"], sep = "")))
    fastalines = c(fastalines, as.character(data[rowNum, "sequences"]))
  }
  fileConn <- file(filename)
  writeLines(fastalines, fileConn)
  close(fileConn)
}

writeFASTA(ABC_sequence_hits, "1_HMMER_sequence_out/ABC_pep_hits.fasta")
writeFASTA(ACR_sequence_hits, "1_HMMER_sequence_out/ACR_pep_hits.fasta")
writeFASTA(SMR_sequence_hits, "1_HMMER_sequence_out/SMR_pep_hits.fasta")
writeFASTA(MatE_sequence_hits, "1_HMMER_sequence_out/MatE_pep_hits.fasta")
writeFASTA(MFS_sequence_hits, "1_HMMER_sequence_out/MFS_pep_hits.fasta")
