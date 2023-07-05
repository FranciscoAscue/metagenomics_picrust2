plot_richness(ps_rar, color = "Procedencia", x = "Procedencia", 
              measures = c("Simpson", "Chao1", "Shannon")) + 
  geom_boxplot(aes(fill = Nivelpeso), alpha=.7) + 
  scale_color_manual(values = 
                       c("#90EE90", "#7CCD7C", "#548B54","#FFD39B","#EEC591","#CDAA7D")) +
  scale_fill_manual(values = 
                      c("#90EE90", "#7CCD7C", "#548B54","#FFD39B","#EEC591","#CDAA7D"))











colores_pcoa <- c("#FFDE8B","#FFA88B","#E9967A","#FA8072")

# PCoA con Bray curtis
plot_ordination(pslog, pcoa_bray, color = "Nivelpeso", shape = "Nivelpeso") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_text(aes(label=Nivelpeso),hjust=0, vjust=1,size = 2)+
  geom_point(size=3)+ 
  scale_fill_manual(values=colores_pcoa) + 
  scale_color_manual(values = colores_pcoa)+
  ggtitle("PCoA Bray Curtis")+
  theme_classic() -> plot_bray
plot_bray





pcoa_plot <- function(physeq_object, distancia = as.character(), 
                      ponderado = FALSE, 
                      colores){
  pcoa_ordination<- ordinate(pslog, method = "PCoA", 
                             distance = as.character(distancia),
                             weighted = FALSE)
  evals <- pcoa_ordination$values$Eigenvalues
  
  # PCoA con Bray curtis
  plot_ordination(pslog, pcoa_ordination, 
                  color = "Nivelpeso", shape = "Nivelpeso") +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    # geom_text(aes(label=sample),hjust=0, vjust=1,size = 2)+
    geom_point(size=3)+ 
    scale_fill_manual(values=colores) + 
    scale_color_manual(values = colores)+
    theme_classic() -> plot_pcoa
  return(plot_pcoa)
}

pcoa_plot(pslog, "bray", ponderado = F, colores_pcoa) +
  ggtitle("PCoA Bray Curtis") -> pcoa_bray
pcoa_bray

#----Ordenamiento Unifrac no ponderado----
pcoa_plot(pslog, "unifrac", ponderado = F, colores_pcoa) +
  ggtitle("PCoA Unifrac no ponderado") -> pcoa_unifrac_uw
pcoa_unifrac_uw


colores <- c("#20DE8B","#CCDE8B","#FFDE8B","#FFA88B","#FF6AD5","#FF6AD5","#C874AA","#AD8CFF","#966BFF","#90CFFF","#2980B9")
plot_abrel <- plot_composition(phylum_abrel, otu.sort = "abundance",
                               x.label = "sample", group_by = "Nivelpeso")

plot_abrel

plot_abrel <- plot_composition(phylum_abrel, otu.sort = "abundance", 
                               x.label = "Samples", group_by = "Nivelpeso")


plot_abrel +
  scale_fill_manual("Phylum", values = colores) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = 8),
        legend.position = "bottom") +
  ggtitle("Relative abundance Phylum")