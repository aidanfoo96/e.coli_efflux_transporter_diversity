# USING THE DENDO FUNCTIONS IN R TO PLOT CLUSTERED HEATMAPS
# HERE ALL TRANSPORTERS ARE PLOTTED AND CLUSTERED


efflux_transporter_tablev3

efflux_transporter_gene_distribution <- data.frame(efflux_transporter_tablev3 %>% 
                                          group_by(target, Phylogroup, superfamily) %>%
                                          summarise(n()) %>%
                                          arrange(target))

transporter_gene_dist_matrix <- as.matrix(efflux_transporter_gene_distribution[, -1])
rownames(transporter_gene_dist_matrix) <- efflux_transporter_gene_distribution$target
transporter.dendro <- as.dendrogram(hclust(d = dist(x = transporter_gene_dist_matrix)))
transporter.dendro.plot <- ggdendrogram(data = transporter.dendro, rotate = TRUE)
transporter.dendro.plot <- transporter.dendro.plot + theme(axis.text.y = element_text(size = 4))

transporter.order <- order.dendrogram(transporter.dendro)
efflux_transporter_gene_distribution$target <- factor(x = efflux_transporter_gene_distribution$target, 
                                           levels = efflux_transporter_gene_distribution$target[transporter.order], 
                                           ordered = TRUE)

transporter_heatmap <- ggplot(data = efflux_transporter_gene_distribution, aes(x = superfamily, y = target)) +
  geom_tile(aes(fill = Phylogroup)) +
  theme(axis.text.y = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "top") 

pdf("MFS_heatmap_clustered.pdf")
grid.newpage()
print(MFS_heatmap, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(MFS.dendro.plot, vp = viewport(x = 0.90, y = 0.465, width = 0.2, height = 1.04))
dev.off()