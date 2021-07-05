# PLOTTING TRUNCATED TRANSPORTERS 

efflux_transporter_tablev3$supergroup <- paste(efflux_transporter_tablev3$superfamily, 
                                               efflux_transporter_tablev3$Phylogroup)

truncation_distribution <- efflux_transporter_tablev3 %>% 
  group_by(supergroup, truncation_assignments) %>%
  summarise(n() / Num_genomes) 


truncation_distribution_plot <- ggplot(truncation_distribution) + 
  aes(x = supergroup, fill = truncation_assignments, y = `n()/Num_genomes`) + 
  geom_bar(position = "fill", stat = "identity") + 
  coord_flip() + 
  labs(x = "", y = "Normalised Gene Transporter Proportion", fill  = "Gene Distribution") + 
  theme_tufte() + 
  scale_fill_manual(values = cbbPalette)

pdf("FIGURES/truncation_distribution.pdf")
print(truncation_distribution_plot)
dev.off()

g <- efflux_transporter_tablev3 %>%
  group_by(Continents) %>%
  summarise(n())
