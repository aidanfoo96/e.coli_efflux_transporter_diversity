plot_heatmap_transporter_proportions <- function(joined_proteins){
  total_num_genes <- joined_proteins %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(total_gene_num = length(gene))
  
  discrete_num_genes <- joined_proteins %>%
    dplyr::group_by(gene, Phylogroup) %>%
    dplyr::summarise(dist_gene_num = n_distinct(gene))
  
  joined_table <- discrete_num_genes %>%
    left_join(total_num_genes, by = "gene") %>%
    mutate(prop_gene_num = dist_gene_num / total_gene_num)
  
  ggplot(joined_table, aes(x = Phylogroup, y = gene)) + 
    geom_tile(aes(fill = prop_gene_num), colour = "white") +
    scale_fill_gradient(low = "blue", high = "green") + 
    theme_classic() 
}

# ABC Transporters 
pdf("efflux_transporter_counts_per_gene_class.pdf")
efflux_transporter_tablev3 %>%
  filter(superfamily == "ATP Binding Cassette (ATP)") %>%
  group_by(target, gene_class) %>%
  summarise(n()) %>%
  rename(count = `n()`) %>%
  arrange(count) %>%
  ggplot() + 
  aes(x = gene_class, y = target) + 
  geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_tufte() + 
  theme(axis.text.y = element_text(size = 5))   
dev.off( )  

# MFS Transporters
pdf("MFS_transporter_counts_per_gene_class.pdf")
efflux_transporter_tablev3 %>%
  filter(superfamily == "Major Facilitator Superfamily (MFS)") %>%
  group_by(target, gene_class) %>%
  summarise(n()) %>%
  rename(count = `n()`) %>%
  arrange(count) %>%
  ggplot() + 
  aes(x = gene_class, y = target) + 
  geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_tufte() + 
  theme(axis.text.y = element_text(size = 5)) 
dev.off()  


# MatE Transporters 
pdf("MatE_transporter_counts_per_gene_class.pdf")
efflux_transporter_tablev3 %>%
  filter(superfamily == "Multi Antimicrobial Extrusion (MATE)") %>%
  group_by(target, gene_class) %>%
  summarise(n()) %>%
  rename(count = `n()`) %>%
  arrange(count) %>%
  ggplot() + 
  aes(x = gene_class, y = target) + 
  geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_tufte() + 
  theme(axis.text.y = element_text(size = 5)) 
dev.off()  

# MDR Transporters 
pdf("MDR_transporter_counts_per_gene_class.pdf")
efflux_transporter_tablev3 %>%
  filter(superfamily == "Multi Drug Resistance (MDR)") %>%
  group_by(target, gene_class) %>%
  summarise(n()) %>%
  rename(count = `n()`) %>%
  arrange(count) %>%
  ggplot() + 
  aes(x = gene_class, y = target) + 
  geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_tufte() + 
  theme(axis.text.y = element_text(size = 5)) 
dev.off()  

# RND Transporters
pdf("RND_transporter_counts_per_gene_class.pdf")
efflux_transporter_tablev3 %>%
  filter(superfamily == "Resistance Nodulation Cell Division (RND)") %>%
  group_by(target, gene_class) %>%
  summarise(n()) %>%
  rename(count = `n()`) %>%
  arrange(count) %>%
  ggplot() + 
  aes(x = gene_class, y = target) + 
  geom_tile(aes(fill = count), colour = "white") + 
  scale_fill_gradient(low = "blue", high = "green") + 
  theme_tufte() + 
  theme(axis.text.y = element_text(size = 5)) 
dev.off() 
