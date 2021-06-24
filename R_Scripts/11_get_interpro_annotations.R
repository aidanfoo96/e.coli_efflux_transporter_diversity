library(tidyverse)
lineage_summary <- read_csv("metadata/F2_lineage_summary.csv")

column_names_interpro <- c("qseqid", "misc", "site", "website", "id", "domain", "two", "three", "evalue", "date", "interproID", "interprodomain", "site2")
ABC_interpro_annotation <- read_tsv("annotation_files/ABC_clustered.fasta.tsv", col_names = column_names_interpro)

column_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
ABC_proteins <- read_tsv("BLAST-BLAST_hits/ABC_clustered_BLAST.tsv", col_names = column_names)

filter_eval <- function(BLAST_table){
  BLAST_table %>%
    filter(evalue < 0.01)
  
}
ABC_proteins <- filter_eval(ABC_proteins)

ABC_proteins_sep <- ABC_proteins %>%
  filter(evalue < 1e-40) %>%
  mutate(gene = qseqid) %>%
  mutate(gene2 = sseqid) %>%
  separate(gene, into = c("id", "gene")) %>%
  separate(gene2, into = c("id2", "gene2")) %>%
  select(qseqid, sseqid, id, id2, evalue) 

ABC_proteins_sep[,3] <- apply(ABC_proteins_sep[, 3], 2, function(x) as.numeric(as.character(x)))
ABC_proteins_sep[,4] <- apply(ABC_proteins_sep[, 4], 2, function(x) as.numeric(as.character(x)))
ABC_interpro_annotation[,9] <- apply(ABC_interpro_annotation[,9], 2, function(x) as.numeric(as.character(x)))
ABC_interpro_annotation$qseqid <- substr(ABC_interpro_annotation$qseqid, 1, nchar(ABC_interpro_annotation$qseqid)-2)
ABC_proteins_sep$qseqid <- substr(ABC_proteins_sep$qseqid, 1, nchar(ABC_proteins_sep$qseqid)-2)
ABC_proteins_sep$sseqid <- substr(ABC_proteins_sep$sseqid, 1, nchar(ABC_proteins_sep$sseqid)-2)


ABC_interpro_annotation2 <- ABC_interpro_annotation %>%
  filter(website == "Pfam") %>%
  filter(domain != "-") %>%
  mutate(sseqid = qseqid) %>%
  select(qseqid, domain)

ABC_interpro_annotation3 <- ABC_interpro_annotation %>%
  filter(website == "Pfam") %>%
  filter(domain != "-") %>%
  mutate(sseqid = qseqid) %>%
  select(sseqid, domain)

ABC_proteins_joiend1 <- ABC_interpro_annotation2 %>%
  left_join(ABC_proteins_sep, by = c("qseqid" = "qseqid")) 

ABC_proteins_joiend2 <- ABC_interpro_annotation3 %>%
  left_join(ABC_proteins_joiend1, by = c("sseqid" = "sseqid")) %>%
  select(qseqid, domain.x, sseqid, domain.y, evalue)

ABC_proteins_joiend3 <- subset(ABC_proteins_joiend2, !duplicated(subset(ABC_proteins_joiend2, select=c(qseqid, sseqid))))

ABC_proteins_joiend2 %>%
  group_by(sseqid, qseqid) %>%
  summarise(domain.y)




ABC_proteins_joiend3 <- lineage_summary %>%
  left_join(ABC_proteins_joiend2, by = c("id" = "id")) %>%
  select(qseqid, sseqid, domain.x, domain.y, Phylogroup, id2)

ABC_proteins_joiend4 <- lineage_summary %>%
  left_join(ABC_proteins_joiend3, by = c("id" = "id2")) %>%
  select(qseqid, sseqid, domain.x, domain.y, Phylogroup.x, Phylogroup.y) 


ABC_proteins_joiend5 <- subset(ABC_proteins_joiend4, !duplicated(subset(ABC_proteins_joiend4, select=c(qseqid, sseqid, domain.y))))

gene_colours <- c("#66FF66", "#FF0066", "#000000", "#FFFF33", "#0000FF", "#009999", 
                  "#00CCFF", "#00CCCC", "#FFCCCC", "#CCFFCC", "#FFFF33", "#00CCFF", 
                  "#666600", "#666600", "#FFFFFF")

names(gene_colours) <- levels(ABC_proteins_joiend5$domain.y)
gene_colours_custom <- scale_colour_manual(name = "Domain", values = gene_colours)

ggplot(ABC_proteins_joiend5) + 
  aes(x = sseqid, ..count.., colour = domain.y, fill = "NA") + 
  geom_bar()  +
  gene_colours_custom                        
  
ggplot(df, aes(Fruit, ..count..)) + geom_bar(aes(fill = Bug), position = "dodge")


write.csv(ABC_proteins_joiend5, 'ABC_proteins_clusterd_with_lineage_and_interpro_annotations2.csv')



filter_eval <- function(BLAST_table){
  BLAST_table %>%
    filter(evalue < 0.01)

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

ABC_network <- convert_for_cytoscape_with_annotv2("BLAST-BLAST_hits/ABC_clustered_BLAST.tsv", "annotation_files/ABC_clustered.fasta.tsv", eval = 1e-50)
MFS_network <- convert_for_cytoscape_with_annotv2("BLAST-BLAST_hits/MFS_clustered_BLAST.csv", "annotation_files/MFS_clustered.fasta.tsv", eval = 1e-100)

write.csv(MFS_network, "MFS_proteins_clustered_with_lineage_and_interpro_annotation.csv")

