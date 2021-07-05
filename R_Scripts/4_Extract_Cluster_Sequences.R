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

cluster_1_60 <- as_tibble(read_csv("ABC_node_list_eval-60_Clust1.csv"))
cluster_1_60_tidy <- tidy_import(cluster_1_60)
ABC_cluster_1_60 <- extract_sequences(cluster_1_60_tidy)
writeFASTAclust(ABC_cluster_1_60, "4_Clustered_Sequences/ABC_cluster1_60.fasta")

cluster_1_50_ng <- read_csv("ABC_node_list_eval-50__Cluster1_nogroup.csv")
cluster_1_50_ng_tidy <- tidy_import(cluster_1_50_ng)
ABC_cluster_1_50_ng <- extract_sequences(cluster_1_50_ng_tidy)
writeFASTAclust(ABC_cluster_1_50_ng, "4_Clustered_Sequences/ABC_cluster1_50_ng.fasta")

cluster_1_40_ng <- read_csv("ABC_node_list_eval-40_nogroup_Cluster1.csv")
cluster_1_40_ng_tidy <- tidy_import(cluster_1_40_ng)
ABC_cluster_1_40_ng <- extract_sequences(cluster_1_40_ng_tidy)
writeFASTAclust(ABC_cluster_1_40_ng, "4_Clustered_Sequences/ABC_cluster1_40_ng.fasta")

cluster_2_40_ng <- read_csv("ABC_node_list_eval-40_Cluster2_nogroup.csv")
cluster_2_40_ng_tidy <- tidy_import(cluster_2_40_ng)
ABC_cluster_2_40_ng <- extract_sequences(cluster_2_40_ng_tidy)
writeFASTAclust(ABC_cluster_2_40_ng, "4_Clustered_Sequences/ABC_cluster2_40_ng.fasta")

cluster_3_40_ng <- read_csv("ABC_node_list_eval-40_Cluster3_nogroup.csv ")
cluster_3_ng_tidy <- tidy_import(cluster_3_40_ng)
ABC_cluster_3_40_ng <- extract_sequences(cluster_3_ng_tidy)
writeFASTAclust(ABC_cluster_3_40_ng, "4_Clustered_Sequences/ABC_cluster3_40_ng.fasta")

cluster_1_55_ng <- read_csv("ABC_node_list_eval-55_Cluster1_nogroup.csv")
cluster_1_55_ng_tidy <- tidy_import(cluster_1_55_ng)
ABC_cluster_1_55_ng <- extract_sequences(cluster_1_55_ng_tidy)
writeFASTAclust(ABC_cluster_1_55_ng, "4_Clustered_Sequences/ABC_cluster1_55_ng.fasta")

cluster_2_55_ng <- read_csv("ABC_node_list_eval-55_nogroup_Cluster2.csv")
cluster_2_55_ng_tidy <- tidy_import(cluster_2_55_ng)
ABC_cluster_2_55_ng <- extract_sequences(cluster_2_55_ng_tidy)
writeFASTAclust(ABC_cluster_2_55_ng, "4_Clustered_Sequences/ABC_cluster2_55_ng.fasta")

cluster1_ABC_efflux <- read_csv("cyto_outs/ACB_clustered_efflux_transporters.csv")
cluster1_ABC_efflux <- tidy_import(cluster1_ABC_efflux)
cluster1_ABC_efflux <- extract_sequences(cluster1_ABC_efflux)
writeFASTAclust(cluster1_ABC_efflux, "4_Clustered_Sequences/cluster1_ABC_efflux.fasta")
writeFASTAclust(cluster1_ABC_efflux, "cyto_outs//cluster1_ABC_efflux.fasta")

cluster2_ABC_efflux <- read_csv("cyto_outs/ABC_transorters_cluster2.csv")
cluster2_ABC_efflux <- tidy_import(cluster2_ABC_efflux)
cluster2_ABC_efflux <- extract_sequences(cluster2_ABC_efflux)
writeFASTAclust(cluster2_ABC_efflux, "4_Clustered_Sequences/cluster2_ABC_efflux.fasta")

cluster3_ABC_efflux <- read_csv("cyto_outs/cluster_3_efflux_transporter_ABC.csv")
cluster3_ABC_efflux <- tidy_import(cluster3_ABC_efflux)
cluster3_ABC_efflux <- extract_sequences(cluster3_ABC_efflux)
writeFASTAclust(cluster3_ABC_efflux, "4_Clustered_Sequences/cluster3_ABC_efflux.fasta")