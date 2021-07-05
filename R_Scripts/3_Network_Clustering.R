library(tidyverse)

# FUNCTION TO FILTER BLAST BLAST HITS TABLE EVALUE TO < 0.01
filter_eval <- function(BLAST_table){
  BLAST_table %>%
    filter(evalue < 0.01)
}

# FUNCTION TO CONVERT BLAST BLAST HITS FOR CYTOSCAPE
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

# FUNCTION TO CONVERT BLAST BLAST HITS FOR CYTOSCAPE WITH ANNOTATIONS FROM INTERPRO
convert_for_cytoscape_with_annot <- function(file, annot_file, eval_threshold){
  column_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  sequence_table <- read_tsv(file, col_names = column_names)
  
  column_names_interpro <- c("qseqid", "misc", "site", "website", "id", "domain", "two", "three", "evalue", "date", "interproID", "interprodomain", "site2")
  interpro_annotation <- read_tsv(annot_file, col_names = column_names_interpro)
  
  sequence_table_filtd <- filter_eval(sequence_table)
  
  sequence_table_filtd$qseqid <- substr(sequence_table_filtd$qseqid, 1, nchar(sequence_table_filtd$qseqid)-2)
  sequence_table_filtd$sseqid <- substr(sequence_table_filtd$sseqid, 1, nchar(sequence_table_filtd$sseqid)-2)
  interpro_annotation$qseqid <- substr(interpro_annotation$qseqid, 1, nchar(interpro_annotation$qseqid)-2)
  
  sequence_table_tidy <- sequence_table_filtd %>%
    select(qseqid, sseqid, evalue, bitscore, pident) %>%
    filter(evalue != 0) %>%
    filter(evalue < eval_threshold)
  
  annotation_tidy <- interpro_annotation %>%
    select(qseqid, website, domain, site2) %>%
    filter(website == c("Pfam", "ProSitePatterns", "SMART", "PANTHER"))
  
  sequences_table_joined <- sequence_table_tidy %>%
    left_join(annotation_tidy, by = "qseqid")
  
  node_list <- sequences_table_joined %>%
    select(qseqid, sseqid, evalue, domain) 
  
  return(node_list)
  
}

ABC_clustered <- convert_for_cytoscape("BLAST-BLAST_hits/run2/ABC_clustered_BLAST.tsv", eval_threshold = 1.00e-50)
MFS_clustered <- convert_for_cytoscape("BLAST-BLAST_hits/run2/MFS_clustered.tsv", eval_threshold = 1.00e-10)
ACR_clustered <- convert_for_cytoscape("BLAST-BLAST_hits/run2/ACR_clustered.tsv", eval_threshold = 1.00e-10)
MatE_clustered <- convert_for_cytoscape("BLAST-BLAST_hits/run2/MatE_clustered.tsv", eval_threshold = 1.00e-3)
SMR_clustered <- convert_for_cytoscape("BLAST-BLAST_hits/run2/SMR_clustered.tsv", eval_threshold = 1.00e-10)


# FUNCTION TO CONVERT BLAST BLAST HITS FOR CYTOSCAPE WITH ANNOTATIONS FROM TCDB .CSV
annotate_with_tcdb <- function(clustered_table){
  protein_family_data <- read_csv("efflux_transporters_tcdb_annotations.csv")
  clustered_table$qseqid2 <- str_sub(clustered_table$qseqid, 3)
  clustered_table$sseqid2 <- str_sub(clustered_table$sseqid, 3)
  clustered_table <- mutate(clustered_table, target = as.character(gsub("^_", "", qseqid2))) 
  clustered_table <- mutate(clustered_table, target2 = as.character(gsub("^_", "", sseqid2))) 
  clustered_table <- clustered_table %>%
    select(qseqid, sseqid, target, target2, evalue)
  
  clustered_table_join1 <- clustered_table %>%
    left_join(protein_family_data, by = c("target" = "target"))
  
  clustered_table_join2 <- clustered_table_join1 %>%
    left_join(protein_family_data, by = c("target2" = "target"))
  
  clustered_table_join2 <- clustered_table_join2 %>%
    select(qseqid, sseqid, evalue, Annotation.x, Annotation.y, Family.x, Family.y)
  
  return(clustered_table_join2)
}

ABC_network_annotated <- annotate_with_tcdb(ABC_clustered)
MFS_network_annotated <- annotate_with_tcdb(MFS_clustered)
ACR_network_annotated <- annotate_with_tcdb(ACR_clustered)
MatE_network_annotated <- annotate_with_tcdb(MatE_clustered)
SMR_network_annotated <- annotate_with_tcdb(SMR_clustered)

# get family information 
write.csv(ABC_network_annotated, 'ABC_node_list_eval-50_nogroup_tcdb_annot.csv')
write.csv(MFS_network_annotated, 'MFS_node_list_eval-10_nogroup_tcdb_annot.csv')
write.csv(ACR_network_annotated, 'ACR_node_list_eval-10_nogroup_tcdb_annot.csv')
write.csv(MatE_network_annotated, 'MatE_node_list_eval-3_nogroup_tcdb_annot.csv')
write.csv(SMR_network_annotated, 'SMR_node_list_eval-10_nogroup_tcdb_annot.csv')


# JUNK : Import files and give column names junk ####
column_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
ABC_sequence_similarities <- read_tsv("BLAST-BLAST_hits/ABC_clustered_BLAST.tsv", col_names = column_names)

column_names_interpro <- c("qseqid", "misc", "site", "website", "id", "domain", "two", "three", "evalue", "date", "interproID", "interprodomain", "site2")
ABC_interpro_annotation <- read_tsv("MFS_clustered.fasta.tsv", col_names = column_names_interpro)
 
MFS_sequence_similarities <- read_tsv("MFS_clustered_BLAST.csv", col_names = column_names)


ABC_interpro_annotation_tidy <- ABC_interpro_annotation %>%
  select(qseqid, website, domain, site2) %>%
  filter(website == c("Pfam", "ProSitePatterns", "SMART", "PANTHER"))

# filter e-values
filter_eval <- function(BLAST_table){
  BLAST_table %>%
    filter(evalue < 0.01)
}


ABC_sequences_filte <- filter_eval(ABC_sequence_similarities)

ABC_sequences_filte$qseqid <- substr(ABC_sequences_filte$qseqid, 1, nchar(ABC_sequences_filte$qseqid)-2)
ABC_sequences_filte$sseqid <- substr(ABC_sequences_filte$sseqid, 1, nchar(ABC_sequences_filte$sseqid)-2)

ABC_sequences_tidy <- ABC_sequences_filte %>%
  select(qseqid, sseqid, evalue, bitscore, pident) %>%
  filter(evalue != 0) %>%
  filter(!grepl("group", qseqid)) %>%
  filter(!grepl("group", sseqid)) %>%
  filter(evalue < 1.00e-40)

ABC_sequences_tidy_joined <- ABC_sequences_tidy %>%
  left_join(ABC_cluster1_tt1, by = "qseqid")

# create vertex and edge list
node_list <- ABC_sequences_tidy %>%
  select(qseqid, sseqid, evalue) 

write.csv(node_list, 'ABC_node_list_eval-40_nogroup.csv')


