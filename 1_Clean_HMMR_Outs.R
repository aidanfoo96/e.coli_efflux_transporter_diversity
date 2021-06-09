# WHAT DOES THIS SCRIPT DO? 
# -- CLEAN THE HMMER TABLES --
# -- GROUP AND PLOT THE HMMER HITS PER SUPERFAMILY -- 


# load the functions scripts before running
# note, make sure the script pathway to this is defined to your wd

rm( list=ls() )
source("scripts/functions_script_ecoli_efflux_diversity_analysis.R")

# Tidy and Group HMMER hits ####
col_names <- c("target", "#" , "accession", "query_name", 
               "E_value", "score", "bias", "E_value2", 
               "score2", "bias2", "exp", "reg", "clu",
               "ov", "env", "dom", "rep", "inc", "description")

MatE_proteins_cut <- as_tibble(read_table2("HMMER_out/MatE_protein_table_ga", col_names = col_names))
ACR_proteins_cut <- as_tibble(read_table2("HMMER_out/ACR_protein_table_ga", col_names = col_names))
ABC_proteins_cut <- as_tibble(read_table2("HMMER_out/ABC_protein_table_ga", col_names = col_names))
MFS_proteins_cut <- as_tibble(read_table2("HMMER_out/MFS_protein_table_ga", col_names = col_names))
EamA_proteins_cut <- as_tibble(read_table2("HMMER_out/EamA_protein_table_ga", col_names = col_names))
MDR_proteins_cut <- as_tibble(read_table2("HMMER_out/MDR_protein_table_ga", col_names = col_names ))

# Gene Frequencies
gene_freqs <- as_tibble(read_csv("metadata/F5_freqs.csv"))

# metadata
ecoli_met <- as_tibble(read_csv("metadata/F1_genome_metadata.csv"))

# lineages summary
ecoli_lin <- as_tibble(read_csv("metadata/F2_lineage_summary.csv"))


# Rename columns 
protein_list <- list(ABC_proteins_cut, ACR_proteins_cut, 
             EamA_proteins_cut, MDR_proteins_cut, 
             MFS_proteins_cut, MatE_proteins_cut)

# Tidy Tables 
protein_list_clean <- lapply(protein_list, clean_the_tablev2)
protein_list_clean2 <- lapply(protein_list, clean_the_table_simple)

# Get protein Tables
ABC_proteins_cut_clean <- protein_list_clean[[1]]
ACR_proteins_cut_clean <- protein_list_clean[[2]]
EamA_proteins_cut_clean <- protein_list_clean[[3]]
MDR_proteins_cut_clean <- protein_list_clean[[4]]
MFS_proteins_cut_clean <- protein_list_clean[[5]]
MatE_proteins_cut_clean <- protein_list_clean[[6]]

ABC_proteins_cut_clean2 <- protein_list_clean2[[1]]
ACR_proteins_cut_clean2 <- protein_list_clean2[[2]]
EamA_proteins_cut_clean2 <- protein_list_clean2[[3]]
MDR_proteins_cut_clean2 <- protein_list_clean2[[4]]
MFS_proteins_cut_clean2 <- protein_list_clean2[[5]]
MatE_proteins_cut_clean2 <- protein_list_clean2[[6]]

# convert chr columns to numeric

i <- 3

ABC_proteins_cut_clean[, i] <- apply(ABC_proteins_cut_clean[, i], 2, function(x) as.numeric(as.character(x)))
ACR_proteins_cut_clean[, i] <- apply(ACR_proteins_cut_clean[, i], 2, function(x) as.numeric(as.character(x)))
EamA_proteins_cut_clean[, i] <- apply(EamA_proteins_cut_clean[, i], 2, function(x) as.numeric(as.character(x)))
MDR_proteins_cut_clean[, i] <- apply(MDR_proteins_cut_clean[, i], 2, function(x) as.numeric(as.character(x)))
MFS_proteins_cut_clean[, i] <- apply(MFS_proteins_cut_clean[, i], 2, function(x) as.numeric(as.character(x)))
MatE_proteins_cut_clean[, i] <- apply(MatE_proteins_cut_clean[, i], 2, function(x) as.numeric(as.character(x)))

# add superfamily
MatE_proteins_superfam <- add_superfamily(MatE_proteins_cut_clean, "Multi Antimicrobial Extrusion (MATE)")
ACR_proteins_superfam <- add_superfamily(ACR_proteins_cut_clean, "Resistance Nodulation Cell Division (RND)")
ABC_proteins_superfam <- add_superfamily(ABC_proteins_cut_clean, "ATP Binding Cassette (ATP)")
MFS_proteins_superfam <- add_superfamily(MFS_proteins_cut_clean, "Major Facilitator Superfamily (MFS)")
EamA_proteins_superfam <- add_superfamily(EamA_proteins_cut_clean, "Small Multidrug Resistance (SMR)")
MDR_proteins_superfam <- add_superfamily(MDR_proteins_cut_clean, "Small Multidrug Resistance (SMR)")

# make one table of transporters
efflux_transporter_superfam <- rbind(MatE_proteins_superfam, 
                                    ACR_proteins_superfam, 
                                    ABC_proteins_superfam,
                                    MFS_proteins_superfam,
                                    EamA_proteins_superfam,
                                    MDR_proteins_superfam)


# summarise the number of lineages a gene is present in 
transporter_gene_dist <- data.frame(efflux_transporter_superfam %>% 
                                          group_by(target, superfamily) %>%
                                          summarise(count = n()))


# steps to make a clustered heat map ####
transporter_gene_dist_matrix <- as.matrix(transporter_gene_dist[, -1])
rownames(transporter_gene_dist_matrix) <- transporter_gene_dist$target
dendro <- as.dendrogram(hclust(d = dist(x = transporter_gene_dist_matrix)))
dendro.plot <- ggdendrogram(data = dendro, rotate = TRUE)
dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 4))

order <- order.dendrogram(dendro)
transporter_gene_dist$target <- factor(x = transporter_gene_dist$target, 
                                           levels = transporter_gene_dist$target[order], 
                                           ordered = TRUE)

# write the heatap to a PDF in your wd
pdf("transporter_heatmap_total.pdf")
ggplot(data = transporter_gene_dist, aes(x = superfamily, y = target)) +
  geom_tile(aes(fill = count)) +
  scale_fill_gradient(low = "#FF7F50", high = "#66CD00") +
  theme(axis.text.x.bottom = element_blank(), 
        legend.position = "left",
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 5)) + 
  coord_flip() + 
  labs(fill = "Lineage Gene Count")
dev.off()

# write a bar chart showing counts of each transporter to wd
pdf("gene_dist_bar.pdf")
efflux_transporter_superfam %>% 
  group_by(superfamily) %>%
  summarise(count = n()) %>%
  ggplot() + 
  aes(x = superfamily, y = count) +
  geom_col(fill = "black") + 
  coord_flip() + 
  theme_tufte() + 
  labs(x = "Gene Count")
dev.off()

# Get tables of the Gene Transporter Distribution
efflux_transporter_superfam2 <- efflux_transporter_superfam %>%
  group_by(superfamily) %>%
  summarise(Gene_Count = n()) %>%
  rename('Superfamily' = "superfamily", 
         GeneCount = "Gene_Count")

pdf("efflux_transporter_table.pdf")
gt(efflux_transporter_superfam2) %>%
  tab_header(title = md("**HMMER Counts Per Superfamily**")) 


# Plot Metadata (extra) #### 

ggplot(ecoli_lin) + 
  aes(x = id, y = Median_num_genes) + 
  geom_col()

ggplot(ecoli_lin) + 
  aes(x = id, y = Num_genomes) + 
  geom_col()

    # very heavy scew 
    # first lineage have majority of genomes
    # consider this bias in further analysis