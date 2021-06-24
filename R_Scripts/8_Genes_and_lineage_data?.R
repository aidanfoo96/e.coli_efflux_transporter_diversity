# Playing the metadata!
library(tidyverse)
library(vcd)
library(scales)
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

clean_the_table <- function(HMMERhits_table){
  HMMERhits_table_filter <- HMMERhits_table %>%
    filter(target != "#") %>%
    filter(target != "#-------------------")
  
  HMMERhits_table_filter$target2 <- HMMERhits_table_filter$target
  HMMERhits_table_filter$target <- str_extract(HMMERhits_table_filter$target, "^.{2}")
  
  HMMERhits_table_filter <- HMMERhits_table_filter %>%
    separate(target2, into = c("id", "gene")) 
  
  HMMERhits_table_filter2 <- HMMERhits_table_filter %>%
    mutate(target = as.character(gsub("_", "", target))) %>%
    select(gene, target, E_value)
  
  HMMERhits_table_filter2[,1] <- apply(HMMERhits_table_filter2[, 1], 2, function(x) as.numeric(as.character(x)))
  
  return(HMMERhits_table_filter2)
  
}

MatE_proteins <- clean_the_table(MatE_proteins_cut)
ACR_proteins <- clean_the_table(ACR_proteins_cut)
ABC_proteins <- clean_the_table(ABC_proteins_cut)
MFS_proteins <- clean_the_table(MFS_proteins_cut)
EamA_proteins <- clean_the_table(EamA_proteins_cut)
MDR_proteins <- clean_the_table(MDR_proteins_cut)

# join lineage summary with genes 
join_protein_hits_with_lineage <- function(tidy_protein_table){
  lineage_summary <- read_csv("metadata/F2_lineage_summary.csv")
  
  protein_lineage_table <- lineage_summary %>%
    inner_join(tidy_protein_table, by = "id")
  
  return(protein_lineage_table)
}

MatE_proteins_lineage_joined <- join_protein_hits_with_lineage(MatE_proteins)
ACR_proteins_lineage_joined <- join_protein_hits_with_lineage(ACR_proteins)
ABC_proteins_lineage_joined <- join_protein_hits_with_lineage(ABC_proteins)
MFS_proteins_lineage_joined <- join_protein_hits_with_lineage(MFS_proteins)
EamA_proteins_lineage_joined <- join_protein_hits_with_lineage(EamA_proteins)
MDR_proteins_lineage_joined <- join_protein_hits_with_lineage(MDR_proteins)

# add superfamily annotation 
add_superfamily <- function(joined_table, superfamily){
  anot <- joined_table %>%
    mutate(superfamily = superfamily)
}

MatE_proteins_lineage_joined2 <- add_superfamily(MatE_proteins_lineage_joined, "Multi Antimicrobial Extrusion (MATE)")
ACR_proteins_lineage_joined2 <- add_superfamily(ACR_proteins_lineage_joined, "Resistance Nodulation Cell Division (RND)")
ABC_proteins_lineage_joined2 <- add_superfamily(ABC_proteins_lineage_joined, "ATP Binding Cassette (ATP)")
MFS_proteins_lineage_joined2 <- add_superfamily(MFS_proteins_lineage_joined, "Major Facilitator Superfamily (MFS)")
EamA_proteins_lineage_joined2 <- add_superfamily(EamA_proteins_lineage_joined, "Resistance Nodulation Cell Division (RND)")
MDR_proteins_lineage_joined2 <- add_superfamily(MDR_proteins_lineage_joined, "Multi Drug Resistance (MDR)")

efflux_transporter_table <- rbind(MatE_proteins_lineage_joined2, 
                                  ACR_proteins_lineage_joined2, 
                                  ABC_proteins_lineage_joined2,
                                  MFS_proteins_lineage_joined2,
                                  EamA_proteins_lineage_joined2,
                                  MDR_proteins_lineage_joined2)




MatE_proteins_lineage_joined %>%
  inner_join(ACR_proteins_lineage_joined)
# which phylogroup? Any patterns? 
plot_stacked_bar_proportions <- function(data, xvar, fillvar){
  data %>%
    ggplot() + 
    aes_string(x = xvar, fill = fillvar) + 
    geom_bar(position = "fill") + 
    labs(y = "Proportion") + 
    theme_minimal()
}

plot_stacked_bar_proportions(MatE_proteins_lineage_joined, 
                                    "Phylogroup", "gene")
plot_stacked_bar_proportions(ACR_proteins_lineage_joined, 
                             "Phylogroup", "gene")
plot_stacked_bar_proportions(ABC_proteins_lineage_joined, 
                             "Phylogroup", "gene")
 

p <- ggplot(ABC_proteins_lineage_joined, 
            aes(x = Phylo))

p <- ABC_proteins_lineage_joined %>%
  group_by(id) %>%
  summarise(mean_genome_number = mean(Num_genomes))


p2 <- ABC_proteins_lineage_joined %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(total_gene_num = length(gene))

p3 <- ABC_proteins_lineage_joined %>%
  dplyr::group_by(gene, Phylogroup) %>%
  dplyr::summarise(dist_gene_num = n_distinct(gene), .groups = "keep")

p4 <- p3 %>%
  left_join(p2, by = "gene") %>%
  mutate(prop_gene_num = dist_gene_num / total_gene_num) 
    # gene proportions of each transporter in each lineage


ggplot(p4, aes(x = Phylogroup, y = gene)) + 
  geom_tile(aes(fill = prop_gene_num), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "red") + 
  theme_classic()

p4 <- p4 %>%
  select(Phylogroup, gene, prop_gene_num)
p5 <- apply(as.matrix(p4), 
            2, 
            as.numeric)

p5 <- data.matrix(p4)


class(p5)

heatmap(p5)

p3 %>% 
  ggplot() + 
  aes(x = gene, y = num_lineages) +
  geom_col() + 
  coord_flip()

data <- read_csv("metadata/F5_freqs.csv")
