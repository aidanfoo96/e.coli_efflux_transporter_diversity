library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
library("ggthemes")
library("tidyverse")

# Make Clustered Heatmap Function ####
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
  
  grid.newpage()
  print(heatmap, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))

}

make_clustered_dendrogram <- function(transporter_gene_list){
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
  
  print(dendro.plot, vp = viewport(x = 0.80, y = 0.465, width = 0.2, height = 1.04))
  
}

ABC_transporters <- efflux_transporter_tablev3 %>%
  filter(superfamily == "ATP Binding Cassette (ATP)") %>%
  group_by(target, gene_class) %>%
  summarise(count = n()) %>%
  arrange(count) 

MatE_transporters <- efflux_transporter_tablev3 %>%
  filter(superfamily == "Multi Antimicrobial Extrusion (MATE)") %>%
  group_by(target, gene_class) %>%
  summarise(count = n()) %>%
  arrange(count) 

MFS_transporters <- efflux_transporter_tablev3 %>%
  filter(superfamily == "Major Facilitator Superfamily (MFS)") %>%
  group_by(target, gene_class) %>%
  summarise(count = n()) %>%
  arrange(count) 

list(MFS_transporters %>% filter(gene_class == "Core" && count < 5))

RND_transporters <- efflux_transporter_tablev3 %>%
  filter(superfamily == "Resistance Nodulation Cell Division (RND)") %>%
  group_by(target, gene_class) %>%
  summarise(count = n()) %>%
  arrange(count) 

SMR_transporters <- efflux_transporter_tablev3 %>%
  filter(superfamily == "Small Multidrug Resistance (SMR)") %>%
  group_by(target, gene_class) %>%
  summarise(count = n()) %>%
  arrange(count) 

ABC_transporters_heat <- make_clustered_heatmap(ABC_transporters)
ABC_transporters_dendo <- make_clustered_dendrogram(ABC_transporters)

SMR_transporters_heat <- make_clustered_heatmap(SMR_transporters)
SMR_transporters_dendo <- make_clustered_dendrogram(SMR_transporters)

MFS_transporters_heat <- make_clustered_heatmap(MFS_transporters)
MFS_transporters_dendo <- make_clustered_dendrogram(MFS_transporters)

RND_transporters_heat <- make_clustered_heatmap(RND_transporters)
RND_transporters_dendo <- make_clustered_dendrogram(RND_transporters)

MATE_transporter_heat <- make_clustered_heatmap(MatE_transporters)
MATE_transporter_dendo <- make_clustered_dendrogram(MatE_transporters)


pdf("FIGURES/ABC_heatmap2.pdf")
print(ABC_transporters_heat)
dev.off()
pdf("FIGURES/ABC_dendo.pdf")
print(ABC_transporters_dendo)
dev.off()

pdf("FIGURES/MatE_Heatmap2.pdf")
print(MATE_transporter_heat)
dev.off()
pdf("FIGURES/MatE_dendo.pdf")
print(MATE_transporter_dendo)
dev.off()

pdf("FIGURES/MFS_heatmap.pdf")
print(MFS_transporters_heat)
dev.off()
pdf("FIGURES/MFS_dendo.pdf")
print(MFS_transporters_dendo)
dev.off()

pdf("FIGURES/RND_heatmap.pdf")
print(RND_transporters_heat)
dev.off()
pdf("FIGURES/RND_dendo.pdf")
print(RND_transporters_dendo)
dev.off()

pdf("FIGURES/SMR_heatmap.pdf")
print(SMR_transporters_heat)
dev.off()

pdf("FIGURES/SMR_dendo.pdf")
print(SMR_transporters_dendo)
dev.off()



# Make Tables showing number of genes in each gene_class ####
library(gt)
count_gene_transporter_class <- function(transporter_table){
  grouped_transporter_table <- transporter_table %>%
    group_by(gene_class) %>%
    summarise(Count = n()) %>%
    rename(GeneClass = "gene_class")
  
  return(gt(grouped_transporter_table))
}

ABC_gene_transporter_class_count <- count_gene_transporter_class(ABC_transporters)
MatE_gene_transporter_class_count <- count_gene_transporter_class(MatE_transporters)
MFS_gene_transporter_class_count <- count_gene_transporter_class(MFS_transporters)
RND_gene_transporter_class_count <- count_gene_transporter_class(RND_transporters)
SMR_gene_transporter_class_count <- count_gene_transporter_class(SMR_transporters)


# stacked bar charts showing number of transporters in each phylogroup 
phylogroup_distribution <- efflux_transporter_tablev3 %>%
  group_by(superfamily, Phylogroup) %>%
  summarise(gene_count = n())

phylogroup_distribution_plot <- ggplot(phylogroup_distribution) + 
  aes(x = Phylogroup, fill = superfamily, y = gene_count) + 
  geom_bar(position = "dodge", stat = "identity") + 
  coord_flip() + 
  labs(x = "", y = "Mean Transporter Count", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)

pdf("phylogroup_distribution.pdf")
print(phylogroup_distribution_plot)
dev.off()
### Junk Script For Reference ####
MFS_transporter_gene_dist <- data.frame(efflux_transporter_tablev3 %>% 
                                          filter(superfamily == "Major Facilitator Superfamily (MFS)") %>%
                                          group_by(target, gene_class) %>%
                                          summarise(n()) %>%
                                          rename(count = `n()`) %>%
                                          arrange(gene_class)) 

MFS_transporter_gene_dist_matrix <- as.matrix(MFS_transporter_gene_dist[, -1])
rownames(MFS_transporter_gene_dist_matrix) <- MFS_transporter_gene_dist$target
MFS.dendro <- as.dendrogram(hclust(d = dist(x = MFS_transporter_gene_dist_matrix)))
MFS.dendro.plot <- ggdendrogram(data = MFS.dendro, rotate = TRUE)
MFS.dendro.plot <- MFS.dendro.plot + theme(axis.text.y = element_text(size = 4))

MFS.order <- order.dendrogram(MFS.dendro)
MFS_transporter_gene_dist$target <- factor(x = MFS_transporter_gene_dist$target, 
                                           levels = MFS_transporter_gene_dist$target[MFS.order], 
                                           ordered = TRUE)
MFS_heatmap <- ggplot(data = MFS_transporter_gene_dist, aes(x = gene_class, y = target)) +
  geom_tile(aes(fill = count)) +
  scale_fill_gradient(low = "#FF7F50", high = "#66CD00") +
  theme(axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "top") 

pdf("MFS_heatmap_clustered.pdf")
grid.newpage()
print(MFS_heatmap, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(MFS.dendro.plot, vp = viewport(x = 0.90, y = 0.465, width = 0.2, height = 1.04))
dev.off()

dfc <- cbind(MFS_transporter_gene_dist, 
             id = seq(nrow(MFS_transporter_gene_dist)),
             cluster = k$cluster)

ggplot(MFS_transporter_gene_dist) + 
  aes(x = gene_class, y = target) + 
  geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_tufte() + 
  theme(axis.text.y = element_text(size = 5)) 