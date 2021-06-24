#---How To USE---#
  # main functions for the majority of the E. coli analysis
  # Aidan Foo 2021 - Liverpool School of Tropical Medicine - 


# Functions
library(tidyverse)
library(Biostrings)
library(ggdendro)
library(reshape)
library(grid)
library(ggthemes)
# Manipulate HMMER outputs ####

# Function to clean HMMER outputs.Produces a tibble with an id column to represent the lineage 
# and a target column to represent the gene hit with lineage variation information retained
clean_the_tablev2 <- function(HMMERhits_table){
  HMMERhits_table_filter <- HMMERhits_table %>%
    filter(target != "#") %>%
    filter(target != "#-------------------")
  
  HMMERhits_table_filter %>%
    select(target, accession, E_value)
  
  HMMERhits_table_filter$target2 <- HMMERhits_table_filter$target
  HMMERhits_table_filter$target <- substr(HMMERhits_table_filter$target, 1, nchar(HMMERhits_table_filter$target)-2)
  HMMERhits_table_filter$target <- str_sub(HMMERhits_table_filter$target, 3)
  HMMERhits_table_filter <- mutate(HMMERhits_table_filter, target = as.character(gsub("^_", "", target))) 
  
  HMMER_table_final <- HMMERhits_table_filter %>%
    separate(target2, into = c("id", "gene")) %>%
    mutate(id = as.numeric(id)) %>%
    select(id, target, E_value)
  
  
  
}

clean_the_table_simple <- function(HMMERhits_table){
  HMMERhits_table_filter <- HMMERhits_table %>%
    filter(target != "#") %>%
    filter(target != "#-------------------") 
    
  HMMERhits_table_filter %>%
    select(target)
  
  
}

# Add superfamily information 
add_superfamily <- function(joined_table, superfamily){
  anot <- joined_table %>%
    mutate(superfamily = superfamily)
}


# Apply MetaData to HMMER Outputs ####
    # Metadata from Gal Horesh E. coli genome repository
join_protein_hits_with_lineage <- function(tidy_protein_table){
  lineage_summary <- read_csv("metadata/F2_lineage_summary.csv")
  
  protein_lineage_table <- lineage_summary %>%
    inner_join(tidy_protein_table, by = "id")
  
  return(protein_lineage_table)
}

join_with_gene_classification_data <- function(protein_table){
  
  protein_table2 <- protein_table %>%
    left_join(gene_classification, by = c("target" = "gene")) %>%
    na.exclude()
}

add_superfamily <- function(joined_table, superfamily){
  anot <- joined_table %>%
    mutate(superfamily = superfamily)
}
# Convert HMMER hits to Sequence information ####
  #These functions add sequences to gene names

# remove redundant characters from HMMER hits
fix_target_name <- function(table){
  as_tibble(substr(table$target, 1, nchar(table$target)-2)) # remove last two characters
}

# join genes with E. coli sequences
join_with_sequences_ecoli <- function(reference_table, protein_table){
  protein_table %>%
    left_join(reference_table, by = c('target' = 'seq_names'))
}

# do all in one go
extract_sequences_from_e_coli <- function(HMMER_hits){
  
  e_coli_sequences <- readAAStringSet("e_coli_peptide.fa") # load from current working directory
  seq_names <- names(e_coli_sequences)
  sequences <- paste(e_coli_sequences)
  e_coli_sequences_df <- as_tibble(data.frame(seq_names, sequences))
  
  HMMER_sequence_hit <- join_with_sequences_ecoli(e_coli_sequences_df, HMMER_hits)
  
}

# Convert the dataframe to a fasta file
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








# Convert BLAST outputs to files ready for cytoscape ####
filter_eval <- function(BLAST_table){
  BLAST_table %>%
    filter(evalue < 0.01)
}
convert_for_cytoscape <- function(file, eval_threshold){
  column_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  sequence_table <- read_tsv(file, col_names = column_names)
  
  sequence_table_filtd <- filter_eval(sequence_table)
  
  sequence_table_filtd$qseqid <- substr(sequence_table_filtd$qseqid, 1, nchar(sequence_table_filtd$qseqid)-2)
  sequence_table_filtd$sseqid <- substr(sequence_table_filtd$sseqid, 1, nchar(sequence_table_filtd$sseqid)-2)
  
  sequence_table_tidy <- sequence_table_filtd %>%
    select(qseqid, sseqid, evalue, bitscore, pident) %>%
    filter(evalue != 0) %>%
    filter(!grepl("group", qseqid)) %>%
    filter(!grepl("group", sseqid)) %>%
    filter(evalue < eval_threshold)
  
  node_list <- sequence_table_tidy %>%
    select(qseqid, sseqid, evalue) 
  
  return(node_list)
  
}
convert_for_cytoscape_with_annotv2 <- function(file, annot_file, eval){
  column_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  sequence_table <- read_tsv(file, col_names = column_names)
  
  column_names_interpro <- c("qseqid", "misc", "site", "website", "id", "domain", "two", "three", "evalue", "date", "interproID", "interprodomain", "site2")
  interpro_annotation <- read_tsv(annot_file, col_names = column_names_interpro)
  
  sequence_table2 <- sequence_table %>%
    distinct(qseqid, sseqid, .keep_all = TRUE) %>%
    mutate(gene = qseqid) %>%
    mutate(gene2 = sseqid) %>%
    separate(gene, into = c("id", "gene")) %>%
    separate(gene2, into = c("id2", "gene2")) %>%
    filter(evalue < eval)
  select(qseqid, sseqid, id, id2, evalue)
  
  sequence_table2[,3] <- apply(sequence_table2[, 3], 2, function(x) as.numeric(as.character(x)))
  sequence_table2[,4] <- apply(sequence_table2[, 4], 2, function(x) as.numeric(as.character(x)))
  interpro_annotation[,9] <- apply(interpro_annotation[,9], 2, function(x) as.numeric(as.character(x)))
  interpro_annotation$qseqid <- substr(interpro_annotation$qseqid, 1, nchar(interpro_annotation$qseqid)-2)
  interpro_annotation$sseqid <- substr(interpro_annotation$sseqid, 1, nchar(interpro_annotation$sseqid)-2)
  sequence_table2$qseqid <- substr(sequence_table2$qseqid, 1, nchar(sequence_table2$qseqid)-2)
  sequence_table2$sseqid <- substr(sequence_table2$sseqid, 1, nchar(sequence_table2$sseqid)-2)
  
  interpro_annotation2 <- interpro_annotation %>%
    filter(website == "Pfam") %>%
    filter(domain != "-") %>%
    mutate(sseqid = qseqid) %>%
    select(qseqid, domain)
  
  interpro_annotation3 <- interpro_annotation %>%
    filter(website == "Pfam") %>%
    filter(domain != "-") %>%
    mutate(sseqid = qseqid) %>%
    select(sseqid, domain)
  
  proteins_joined <- sequence_table2 %>%
    left_join(interpro_annotation2, by = c("qseqid" = "qseqid")) 
  
  proteins_joined2 <- interpro_annotation3 %>%
    left_join(proteins_joined, by = c("sseqid" = "sseqid")) 
  
  proteins_joined3 <- subset(proteins_joined2, !duplicated(subset(proteins_joined2, select=c(qseqid, sseqid, domain.y))))
} 
  
# Convert Cytoscape Session Clusters into FASTA Files ####
  # Set of functions to conver the cytoscape CSV export into fasta files

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
# Plot Phylogroup - Gene Distribution Heatmaps ####

plot_heatmap_transporter_proportions <- function(joined_proteins){
  total_num_genes <- joined_proteins %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(total_gene_num = length(gene))
  
  discrete_num_genes <- joined_proteins %>%
    dplyr::group_by(gene, Phylogroup) %>%
    dplyr::summarise(dist_gene_num = n_distinct(gene))
  
  joined_table <- discrete_num_genes %>%
    left_join(total_num_genes, by = "gene") %>%
    mutate(prop_gene_num = dist_gene_num / total_gene_num)
  
  ggplot(joined_table, aes(x = Phylogroup, y = gene)) + 
    geom_tile(aes(fill = prop_gene_num), colour = "white") +
    scale_fill_gradient(low = "blue", high = "green") + 
    theme_classic() 
  
}
# Make clustered heatmaps for gene 'rarity' metadata ####
    # insert a gene transporter list in tibble format

make_clustered_heatmap <- function(transporter_gene_list){
  transporter_gene_list_matrix <- as.matrix(transporter_gene_list[, -1])
  rownames(transporter_gene_list_matrix) <- transporter_gene_list$target
  dendro <- as.dendrogram(hclust(d = dist(x = transporter_gene_list_matrix)))
  dendro.plot <- ggdendrogram(data = dendro, rotate = TRUE)
  dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 4))
  
  order <- order.dendrogram(dendro)
  transporter_gene_list$target <- factor(x = transporter_gene_list$target, 
                                         levels = transporter_gene_list$target[order], 
                                         ordered = TRUE)
  heatmap <- ggplot(data = transporter_gene_list, aes(x = gene_class, y = target)) +
    geom_tile(aes(fill = count)) +
    scale_fill_gradient(low = "#FF7F50", high = "#66CD00") +
    theme_tufte() + 
    theme(axis.text.y = element_blank(), 
          axis.title = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = "top")
  
  return(heatmap)
  
}




