library(treeio)
library(ggtree)
library(tidyverse)
library(tidytree)
library(ape)
library(ggstance)

# Load Data ####
MatE_tree <- read.mrbayes("3_MrBayes_Tres/MatE_trim_tree.out.con.tre")
ACR_tree <- read.mrbayes("3_MrBayes_Tres/ACR_pep_tree.out.con.tre")
MDR_tree <- read.mrbayes("3_MrBayes_Tres/MDR_pep_tree.out.con.tre")
EamA_tree <- read.mrbayes("3_MrBayes_Tres/EamA_pep_tree.out.con.tre")
ABC_tree <- read.mrbayes("3_MrBayes_Tres/ABC_cluster1_tree.out.con.tre")
MFS_tree <- read.mrbayes("3_MrBayes_Tres/MFS_trim_tree.out.con.tre")

# Make Tree Daylight Tree ####
make_daylight_tree <- function(mrbayes_tree){
  p <- ggtree(mrbayes_tree, layout = "daylight")
  p$data <- p$data %>%
    separate(label, into = c("lineage", "label"), sep = "_")
  p$data[,4] <- sapply(p$data[,4], as.numeric)
  
  ecoli_lin <- as_tibble(read_csv("metadata/F2_lineage_summary.csv"))
  ecoli_lin2 <- ecoli_lin %>%
    rename(lineage = id) %>%
    left_join(p$data, by = "lineage")
  
  p2 <- p %<+% ecoli_lin2
  
  p2 + geom_tippoint(aes(color = label))
  
}

ACR_tree_dl <- make_daylight_tree(ACR_tree) # RND family
MatE_tree_dl <- make_daylight_tree(MatE_tree)
MDR_tree_dl <- make_daylight_tree(MDR_tree) # SMR family? 
EamA_tree_dl<- make_daylight_tree(EamA_tree) # 
ABC_tree_dl <- make_daylight_tree(ABC_tree)
MFS_tree_dl <- make_daylight_tree(MFS_tree)

make_circular_tree <- function(mrbayes_tree){
  p <- ggtree(mrbayes_tree, layout = "circular")
  p$data <- p$data %>%
    separate(label, into = c("lineage", "label"), sep = "_")
  p$data[,4] <- sapply(p$data[,4], as.numeric)
  
  ecoli_lin <- as_tibble(read_csv("metadata/F2_lineage_summary.csv"))
  ecoli_lin2 <- ecoli_lin %>%
    rename(lineage = id) %>%
    left_join(p$data, by = "lineage")
  
  p2 <- p %<+% ecoli_lin2
  
  p2 + geom_tippoint(aes(color = label))
  
}

ACR_tree_c <- make_circular_tree(ACR_tree) # RND family
MatE_tree_c <- make_circular_tree(MatE_tree)
MDR_tree_c <- make_circular_tree(MDR_tree) # SMR family? 
EamA_tree_c <- make_circular_tree(EamA_tree) # 
ABC_tree_c <- make_circular_tree(ABC_tree)

# Pivot Gene Frequencies #####
get_gene_frequencies <- function(phylogeny){
  gene_freqs <- as_tibble(read_csv("metadata/F5_freqs.csv"))
  q <- ggtree(phylogeny)
  q$data$label <- substr(q$data$label, 1, nchar(q$data$label)-2)
  
  gene_freqs_piv <- gene_freqs %>%
    pivot_longer(cols = `1`:`51`, 
                 names_to = "lineage", 
                 values_to = "freq") %>%
    rename(label = Gene) 
  
  gene_freqs_piv_tidy <- gene_freqs_piv %>%
    unite(gene_freqs_piv, lineage, label, sep = "_") %>%
    rename(label = gene_freqs_piv)
  
  gene_freqs_piv_filt <- gene_freqs_piv_tidy %>%
    filter(label %in% q$data$label) %>%
    separate(label, into = c("lineage", "label"), sep = "_")
  
  return(gene_freqs_piv_filt)
  
}

ACR_gene_freq <- get_gene_frequencies(ACR_tree)
MatE_gene_freq <- get_gene_frequencies(MatE_tree)
MDR_gene_freq <- get_gene_frequencies(MDR_tree)
EamA_gene_freq <- get_gene_frequencies(EamA_tree)
ABC_gene_freq <- get_gene_frequencies(ABC_tree)

plot_gene_frequencies <- function(gene_freqs){
  ggplot(gene_freqs) + 
    aes(x = lineage, y = freq, col = label) + 
    geom_point(size = 0.9) +
    coord_flip() + 
    theme_minimal()
}

plot_gene_frequencies(ACR_gene_freq)
plot_gene_frequencies(MatE_gene_freq)
plot_gene_frequencies(EamA_gene_freq)
plot_gene_frequencies(MDR_gene_freq)
plot_gene_frequencies(ABC_gene_freq)


m <- ggtree(ACR_tree)
m$data$label
# Convert back into Phylogenetic Tree
p <- ggtree(ACR_tree, layout = "daylight")

p$data <- p$data %>%
  separate(label, into = c("lineage", "label"), sep = "_")

p$data[, 4] <- sapply(p$data[, 4], as.numeric)

ecoli_lin2 <- ecoli_lin %>%
  rename(lineage = id) %>%
  left_join(p$data, by = "lineage")

gene_freqs2 <- gene_freqs_piv_tidy

p2 <- p %<+% ecoli_lin2
p2$data

p2 + geom_tippoint(aes(color = label))

ecoli_lin + 
  geom_tiplab()

tree2 <- groupClade(tree, c(1,8))

# Load Lineage Data
ecoli_lin <- as_tibble(read_csv("F2_lineage_summary.csv"))
gene_freqs <- as_tibble(read_csv("F5_freqs.csv"))

# Extract label - use for 
