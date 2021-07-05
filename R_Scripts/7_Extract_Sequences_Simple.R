library(tidyverse)
library(Biostrings)

# --- Input a cytoscape CSV file --- # 
convert_cyto_out_to_fasta <- function(cyto_csv) {
  
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
  
  cluster_file <- as_tibble(read_csv(cyto_csv))
  cluster_file_tidy <- tidy_import(cluster_file)
  cluster_file_tidy_sequences <- extract_sequences(cluster_file_tidy)
  
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


# DO THIS FOR DESIRED CYTOSCAPE FILE
MFS_10_clust1.1 <- convert_cyto_out_to_fasta("cluster_csvs//MFS_node_list_eval-10_nogroup_Cluster1.1.csv")  
writeFASTAclust(MFS_10_clust1.1, "2_clustered_sequences (USEARCH)/Unaligned/MFS_cluster1_10_ng_unaligned.fasta")

ABC_40_clust <- convert_cyto_out_to_fasta("cluster_csvs/ABC_node_list_eval-40_nogroup_Cluster1.csv")
writeFASTAclust(ABC_40_clust, "4_Clustered_Sequences/ABC_Cluster1_ng_unaligned.fasta")

