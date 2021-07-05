plot_heatmap_transporter_proportions <- function(joined_proteins, category){
  total_num_genes <- joined_proteins %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(total_gene_num = length(gene))

  discrete_num_genes <- joined_proteins %>%
  dplyr::group_by(gene, category) %>%
  dplyr::summarise(dist_gene_num = n_distinct(gene))

  joined_table <- discrete_num_genes %>%
  left_join(total_num_genes, by = "gene") %>%
  mutate(prop_gene_num = dist_gene_num / total_gene_num)

  ggplot(joined_table, aes(x = category, y = gene)) + 
  geom_tile(aes(fill = prop_gene_num), colour = "white") +
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_classic() 

}

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
