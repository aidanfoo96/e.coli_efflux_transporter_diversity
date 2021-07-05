library(tidyverse)
library(ggthemes)

gene_classification <- read_tsv("metadata/classification_v2.csv")
col_names <- c("target", "#" , "accession", "query_name", 
               "E_value", "score", "bias", "E_value2", 
               "score2", "bias2", "exp", "reg", "clu",
               "ov", "env", "dom", "rep", "inc", "description")

MatE_proteins_cut <- as_tibble(read_table2("HMMER_out/MatE_protein_table_ga", col_names = col_names))
ACR_proteins_cut <- as_tibble(read_table2("HMMER_out/ACR_protein_table_ga", col_names = col_names))
ABC_proteins_cut <- as_tibble(read_table2("HMMER_out/ABC_protein_table_ga", col_names = col_names))
ABC_membrane <- as_tibble(read_table2("HMMER_out/ABC_membrane_table_ga", col_names = col_names))
ABC_membrane_2 <- as_tibble(read_table2("HMMER_out/ABC_membrane_2_table_ga", col_names = col_names))
ABC2_proteins_cut <- as_tibble(read_table2("HMMER_out/ABC2_membrane_table_ga", col_names = col_names))
ABC2_membrane_3 <- as_tibble(read_table2("HMMER_out/ABC2_membrane_3_table_ga", col_names = col_names))
MFS_proteins_cut <- as_tibble(read_table2("HMMER_out/MFS_protein_table_ga", col_names = col_names))
EamA_proteins_cut <- as_tibble(read_table2("HMMER_out/EamA_protein_table_ga", col_names = col_names))
MDR_proteins_cut <- as_tibble(read_table2("HMMER_out/MDR_protein_table_ga", col_names = col_names ))

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
    select(id, target)
  
  
  
}

protein_list <- list(ABC_proteins_cut, ACR_proteins_cut, 
                     EamA_proteins_cut, MDR_proteins_cut, 
                     MFS_proteins_cut, MatE_proteins_cut, ABC_membrane, ABC_membrane_2, 
                     ABC2_proteins_cut, ABC2_membrane_3)

# Tidy Tables 
protein_list_clean <- lapply(protein_list, clean_the_tablev2)

# Get protein Tables
ABC_proteins_cut_clean <- protein_list_clean[[1]]
ACR_proteins_cut_clean <- protein_list_clean[[2]]
EamA_proteins_cut_clean <- protein_list_clean[[3]]
MDR_proteins_cut_clean <- protein_list_clean[[4]]
MFS_proteins_cut_clean <- protein_list_clean[[5]]
MatE_proteins_cut_clean <- protein_list_clean[[6]]
ABC_membrane_clean <- protein_list_clean[[7]]
ABC_membrane_2_clean <- protein_list_clean[[8]]
ABC2_proteins_cut_clean <- protein_list_clean[[9]]
ABC2_membrane_3_clean <- protein_list_clean[[10]]

# join all the ABC transporters and aggregate to remove duplicate hits
ABC_proteins_full_clean2 <- rbind(ABC_proteins_cut_clean, 
                                  ABC_membrane_clean,
                                  ABC_membrane_2_clean,
                                  ABC2_proteins_cut_clean,
                                  ABC2_membrane_3_clean)

ABC_proteins_full_clean2 <- ABC_proteins_full_clean2 %>%
  distinct(target, id, .keep_all = TRUE)

# join MDR and EamA transporters, aggregate
SMR_transporters_full_clean2 <- rbind(MDR_proteins_cut_clean, 
                                      EamA_proteins_cut_clean)

SMR_proteins_full_clean2 <- SMR_transporters_full_clean2 %>%
  distinct(target, id, .keep_all = TRUE) # from 800 --> 751 hits



add_superfamily <- function(joined_table, superfamily){
  anot <- joined_table %>%
    mutate(superfamily = superfamily)
}

MatE_proteins <- add_superfamily(MatE_proteins, "Multi Antimicrobial Extrusion (MATE)")
ACR_proteins <- add_superfamily(ACR_proteins, "Resistance Nodulation Cell Division (RND)")
ABC_proteins <- add_superfamily(ABC_proteins_full_clean2, "ATP Binding Cassette (ATP)")
MFS_proteins <- add_superfamily(MFS_proteins, "Major Facilitator Superfamily (MFS)")
SMR_proteins <- add_superfamily(SMR_proteins_full_clean2, "Small Multidrug Resistance (SMR)")

join_with_gene_classification_data <- function(protein_table){
  
  protein_table2 <- protein_table %>%
    left_join(gene_classification, by = c("target" = "gene")) %>%
    na.exclude()
}

MatE_proteins_gene_classif <- join_with_gene_classification_data(MatE_proteins)
ACR_proteins_gene_classif <- join_with_gene_classification_data(ACR_proteins)
ABC_proteins_gene_classif <- join_with_gene_classification_data(ABC_proteins)
MFS_proteins_gene_classif <- join_with_gene_classification_data(MFS_proteins)
EamA_proteins_gene_classif <- join_with_gene_classification_data(EamA_proteins)
MDR_proteins_gene_classif <- join_with_gene_classification_data(MDR_proteins)

join_protein_hits_with_lineage <- function(tidy_protein_table){
  lineage_summary <- read_csv("metadata/F2_lineage_summary.csv")
  
  protein_lineage_table <- tidy_protein_table %>%
    inner_join(lineage_summary, by = "id")
  
  return(protein_lineage_table)
}

MatE_proteins_gene_class_lin <- join_protein_hits_with_lineage(MatE_proteins_gene_classif)
ACR_proteins_gene_class_lin <- join_protein_hits_with_lineage(ACR_proteins_gene_classif)
ABC_proteins_gene_class_lin <- join_protein_hits_with_lineage(ABC_proteins_gene_classif)
MFS_proteins_gene_class_lin <- join_protein_hits_with_lineage(MFS_proteins_gene_classif)
EamA_proteins_gene_class_lin <- join_protein_hits_with_lineage(EamA_proteins_gene_classif)
MDR_proteins_gene_class_lin <- join_protein_hits_with_lineage(MDR_proteins_gene_classif)

efflux_transporter_tablev2 <- bind_rows(MatE_proteins_gene_classif, 
                                    ACR_proteins_gene_classif, 
                                    ABC_proteins_gene_classif,
                                    MFS_proteins_gene_classif,
                                    EamA_proteins_gene_classif,
                                    MDR_proteins_gene_classif)


efflux_transporter_tablev3 <- bind_rows(MatE_proteins_gene_class_lin, 
                                    ACR_proteins_gene_class_lin, 
                                    ABC_proteins_gene_class_lin, 
                                    MFS_proteins_gene_class_lin, 
                                    EamA_proteins_gene_class_lin, 
                                    MDR_proteins_gene_class_lin)



## Summary Statistics ####

lineage_gene_statistics <- ABC_proteins_gene_classif %>%
    group_by(id) %>%
    summarise(
      Core = n_distinct(core),
      Intermediate = n_distinct(inter),
      Rare = n_distinct(rare)) %>%
    pivot_longer(!id, names_to = "stat", values_to = "value") 

lin_plot <- ggplot(lineage_gene_statistics) + 
  aes(x = id, fill = stat, y = value) + 
  geom_bar(position = "stack", stat = "identity") + 
  coord_flip() + 
  labs(x = "Lineage", y = "No. ABC Transporter Genes", fill = "Gene Distribution") + 
  theme_tufte()

transporter_gene_distribution <- efflux_transporter_tablev2 %>%
  group_by(superfamily) %>%
  summarise(
    Core = mean(core), 
    Intermediate = mean(inter), 
    Rare = mean(rare)) %>%
  pivot_longer(!superfamily, names_to = "stat", values_to = "value")

pdf("mean_transporter_count_per_superfamily.pdf")
ggplot(transporter_gene_distribution) + 
  aes(x = superfamily, fill = stat, y = value) + 
  geom_bar(position = "dodge", stat = "identity") + 
  coord_flip() + 
  labs(x = "", y = "Mean Transporter Count", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)
dev.off()

transporter_lineage_distribution <- efflux_transporter_tablev2 %>%
  group_by(id) %>%
  summarise(
    mean_core_gene = mean(core), 
    mean_inter_gene = mean(inter), 
    mean_rare_genes = mean(rare)) %>%
  pivot_longer(!superfamily, names_to = "stat", values_to = "value")

efflux_transporter_tablev3$supergroup <- paste(efflux_transporter_tablev3$superfamily, 
                                               efflux_transporter_tablev3$Phylogroup)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

transporter_pathosupergroup_distribution <- efflux_transporter_tablev3 %>%
  group_by(supergroup, gene_class) %>%
  summarise(count = n()) 

transporter_pathosupergroup_distribution2 <- transporter_pathosupergroup_distribution %>%
  left_join(efflux_transporter_tablev3) %>%
  mutate(normalised_counts = count / Num_genomes) %>%
  group_by(supergroup, gene_class) %>%
  summarise(MeanNormalisedCounts = mean(normalised_counts))

pdf("normalised_mean_phylogroup_gene_distribution6.pdf")
ggplot(transporter_pathosupergroup_distribution2) + 
  aes(x = supergroup, fill = gene_class, y = MeanNormalisedCounts) + 
  geom_bar(position = "stack", stat = "identity") + 
  coord_flip() + 
  labs(x = "", y = "Gene Transporter Proportion", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)
dev.off()

pdf("phylogroup_gene_distribution5.pdf")
ggplot(transporter_pathosupergroup_distribution) + 
  aes(x = supergroup, fill = gene_class, y = `n()`) + 
  geom_bar(position = "dodge", stat = "identity") + 
  coord_flip() + 
  labs(x = "", y = "Gene Transporter Proportion", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)
dev.off()


transporter_gene_class_distribution <- efflux_transporter_tablev3 %>%
  group_by(target, gene_class) %>%
  summarise(
    count = n(), 
  ) 

transporter_gene_class_distribution2 <- transporter_gene_class_distribution %>%
  left_join(efflux_transporter_tablev3) %>%
  mutate(normalised_counts = count / Num_genomes) %>%
  group_by(target, gene_class) %>%
  summarise(MeanNormalisedCounts = sum(normalised_counts))

pdf("normalised_transporter_gene_class_count.pdf")
ggplot(transporter_gene_class_distribution2) + 
  aes(x = gene_class, fill = gene_class, y = MeanNormalisedCounts) + 
  geom_bar(position = "stack", stat = "identity") + 
  labs(x = "", y = "Transporter Count", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)
dev.off()

pdf("transporter_gene_class_count.pdf")
ggplot(transporter_gene_class_distribution) + 
  aes(x = gene_class, fill = gene_class, y = count) + 
  geom_bar(position = "stack", stat = "identity") + 
  labs(x = "", y = "Transporter Count", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)
dev.off()


# JUNK NEEDS FURTHER WORK #####

transporter_pathogroup_distribution <- efflux_transporter_tablev3 %>%
  group_by(superfamily, Phylogroup) %>%
  summarise(
    Core = n_distinct(core), 
    Intermediate = n_distinct(inter), 
    Rare = n_distinct(rare)
  ) %>%
  pivot_longer(!c(superfamily, Phylogroup), names_to = "stat", values_to = "value")

transporter_pathogroup_distribution$supergroup <- paste(transporter_pathogroup_distribution$superfamily, 
                                                        transporter_pathogroup_distribution$Phylogroup)

transporter_pathogroup_distributionv2 <- transporter_pathogroup_distribution %>%
  group_by(supergroup) %>%
  summarise(
    Core = n_distinct(core), 
    Intermediate = n_distinct(inter), 
    Rare = n_distinct(rare)
  ) %>%
  pivot_longer(!c(superfamily, Phylogroup), names_to = "stat", values_to = "value")


transporter_pathogroup_distribution <- arrange(transporter_pathogroup_distribution, 
                                               Phylogroup)
transporter_pathogroup_distribution <- as.data.frame(transporter_pathogroup_distribution)

label_data <- transporter_pathogroup_distribution 
  
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$Phylogroup) / number_of_bar

circle_plot <- ggplot(transporter_pathogroup_distribution, 
                      aes(x = as.factor(supergroup), fill = stat, y = value)) + 
  geom_bar(stat = "identity") + 
  ylim(-100,120) + 
  theme_minimal() + 
  theme(
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    panel.grid = element_blank(), 
    plot.margin = unit(rep(-2,4), "cm")
  ) + 
  coord_polar(start = 0)

empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar * nlevels(transporter_pathogroup_distribution$Phylogroup), ncol(transporter_pathogroup_distribution)))
colnames(to_add) <- colnames(transporter_pathogroup_distribution)
to_add$Phylogroup <- rep(levels(transporter_pathogroup_distribution$Phylogroup), each=empty_bar)
transporter_pathogroup_distribution <- rbind(transporter_pathogroup_distribution, to_add)
transporter_pathogroup_distribution <- transporter_pathogroup_distribution %>% arrange(Phylogroup)
transporter_pathogroup_distribution$id <- seq(1, nrow(transporter_pathogroup_distribution))

label_data <- transporter_pathogroup_distribution
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

base_data <- transporter_pathogroup_distribution %>%
  group_by(Phylogroup) %>%
  summarise(start = min(id), end = max(id)) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end)))
  
grid_data <- base_data 
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 0.5
grid_data <- grid_data[-1,]

pdf('phylogroup_plot.pdf')
ggplot(transporter_pathogroup_distribution, aes(x=as.factor(supergroup), y=value, fill=stat)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(supergroup), y=value, fill=stat), stat="identity", position = "fill") +

  geom_bar(aes(x=as.factor(supergroup), y=value, fill=stat), stat="identity", alpha=0.5) +
  ylim(-100,50) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label="", hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -10, xend = end, yend = -10), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -35, label=Phylogroup), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
dev.off()



transporter_pathogroup_distribution <- transporter_pathogroup_distribution %>%
  arrange(Phylogroup, stat)

ggplot(transporter_pathogroup_distribution, aes(x=as.factor(id), y=value, fill=stat)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", position = "dodge", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label="", hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p





