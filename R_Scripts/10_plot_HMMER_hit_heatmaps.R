# heatmaps show the number of transporters within each E. coli lineage

# plot heatmap of each transporter
plot_heatmap_transporter_proportions <- function(joined_proteins){
  total_num_genes <- joined_proteins %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(total_gene_num = length(gene))
  
  discrete_num_genes <- joined_proteins %>%
    dplyr::group_by(gene, Phylogroup) %>%
    dplyr::summarise(dist_gene_num = n_distinct(gene))
  
      # normalisation here
  joined_table <- discrete_num_genes %>%
    left_join(total_num_genes, by = "gene") %>%
    mutate(prop_gene_num = dist_gene_num / total_gene_num)
  
  ggplot(joined_table, aes(x = Phylogroup, y = gene)) + 
    geom_tile(aes(fill = prop_gene_num), colour = "white") +
    scale_fill_gradient(low = "blue", high = "green") + 
    theme_classic() 
  
}

plot_heatmap_transporter_proportions(ABC_proteins_lineage_joined)
plot_heatmap_transporter_proportions(ACR_proteins_lineage_joined)
plot_heatmap_transporter_proportions(MatE_proteins_lineage_joined)
plot_heatmap_transporter_proportions(MFS_proteins_lineage_joined)
plot_heatmap_transporter_proportions(EamA_proteins_lineage_joined)
plot_heatmap_transporter_proportions(MDR_proteins_lineage_joined)

# plot heatmap of all transporters
plot_heatmap_all_transporter_proportions <- function(joined_proteins){
  total_num_genes <- joined_proteins %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(total_gene_num = length(gene))
  
  discrete_num_genes <- joined_proteins %>%
    dplyr::group_by(gene, superfamily, Phylogroup) %>%
    dplyr::summarise(dist_gene_num = n_distinct(gene))
  
  joined_table <- discrete_num_genes %>%
    left_join(total_num_genes, by = "gene") %>%
    mutate(prop_gene_num = dist_gene_num / total_gene_num)
  
  ggplot(joined_table, aes(x = Phylogroup, y = gene)) + 
    geom_tile(aes(fill = superfamily), colour = "white") + 
    scale_fill_manual(values = cbbPalette)
  
}


pdf("efflux_transporter_table.pdf")
plot_heatmap_all_transporter_proportions(efflux_transporter_table)
dev.off()
