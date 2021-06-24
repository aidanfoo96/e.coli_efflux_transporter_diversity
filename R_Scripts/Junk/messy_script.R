library(ggplot2)
library(ggdendro)
library(reshape2)
library(forcats)

total_num_genes <- joined_proteins %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(total_gene_num = length(gene))

discrete_num_genes <- joined_proteins %>%
  dplyr::group_by(gene, superfamily, Phylogroup) %>%
  dplyr::summarise(dist_gene_num = n_distinct(gene))

joined_table <- discrete_num_genes %>%
  left_join(total_num_genes, by = "gene") %>%
  mutate(prop_gene_num = dist_gene_num / total_gene_num)

efflux_transporter_total_genes <- efflux_transporter_table %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(total_gene_num = length(gene))


efflux_transporter_discrete_genes <- efflux_transporter_table %>%
  dplyr::group_by(gene, superfamily, Phylogroup) %>%
  dplyr::summarise(dist_gene_num = length(gene))

efflux_transporters_joined <- efflux_transporter_discrete_genes %>%
  left_join(efflux_transporter_total_genes, by = "gene") %>%
  mutate(prop_gene_num = dist_gene_num / total_gene_num) 
  
efflux_transporters_joined$gene <- with(efflux_transporters_joined, reorder(gene))
  
ggplot(efflux_transporters_joined, aes(x = Phylogroup, y = gene)) + 
  geom_tile(aes(fill = superfamily), colour = "white") 


efflux_trans_dendo <- hclust(d = dist(x = efflux_transporters_matrix))$labels



efflux_transporters_joined$gene <- factor(efflux_transporters_joined$gene, levels = rownames(efflux_trans_dendo)[or])



efflux_transporters_values <- efflux_transporters_joined %>%
  ungroup() %>%
  select(prop_gene_num)

efflux_transporters_matrix <- as.matrix(efflux_transporters_values)
rownames(efflux_transporters_matrix) <- efflux_transporters_joined$gene


efflux_trans_dendo_plot <- ggdendrogram(data = efflux_trans_dendo, rotate = TRUE)

grid.newpage()

transporter_order <- order.dendrogram(efflux_trans_dendo)

efflux_transporters_joined$gene <- factor(x = efflux_transporters_joined$gene, 
                                          levels = efflux_transporters_joined$gene[transporter_order], 
                                          ordered = TRUE)
