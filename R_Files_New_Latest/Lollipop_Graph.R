# Bubble Chart
enrichment_kegg <- read.csv(file = "C:\\Users\\HP\\Documents\\R_Go\\Raw Combined + Filtered\\MT-Genes Filtered Combined\\Latest_1.01_Low_Resolution_C15\\_Downstream_Analysis\\Gene Ontology Report\\GO_Threshold_0.01\\Biological Process/enrichment.csv")
labels <- as.numeric(format(x = seq((min(enrichment_kegg$Enrichment.FDR) - 1e-8), (max(enrichment_kegg$FDR) + 1e-4), length.out = 10), scientific = TRUE, digits = 3))

ggplot(data = enrichment_kegg, mapping = aes(x = `Fold.Enrichment`, y = Term, colour = FDR, size = Count)) + 
    scale_x_continuous(limits = c(0, 5)) + 
    geom_point(alpha = 0.5) + 
    scale_size_continuous(range = c(min(5), max(15))) + 
    # geom_segment(mapping = aes(x = 0, y = Term, xend = `Fold Enrichment`, yend = Term), size = 1) + 
    # geom_point(mapping = aes(size = Count)) + 
    scale_color_gradient(low = "gold", high = "red", limits = c(min(labels), max(labels)), 
                         breaks = labels,
                         labels = as.character(labels),
                         guide = guide_colorbar(
                             title = "FDR", 
                             title.position = "top", 
                             title.hjust = 0.5, 
                             barwidth = unit(1, "lines"), 
                             barheight = unit(15, "lines"), 
                             ticks = TRUE,
                             ticks.colour = "black"
                         ),
                         oob = scales::squish) + 
    theme_bw() + 
    theme(
        axis.text.x = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(t = 15)),
        axis.title.y = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(r = 15)),
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 1),
        legend.title = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(b = 15)),
        legend.text = element_text(size = 15, family = "serif", face = "bold", colour = "black"),
        plot.title = element_text(size = 30, family = "serif", face = "bold", colour = "black"),
        panel.border = element_rect(linewidth = 1)
    ) + 
    labs(
        title = "KEGG Pathway Gene Enrichment Analysis",
        x = "Fold Enrichment",
        y = "Significantly Enriched KEGG Pathways"
    )






# Lollipop Chart
enrichment_kegg <- read.csv(file = "C:\\Users\\HP\\Documents\\R_Go\\Raw Combined + Filtered\\MT-Genes Filtered Combined\\Latest_1.01_Low_Resolution_C15\\_Downstream_Analysis\\DAVID KEGG Pathway Analysis/All_Columns_DAVIDChartReport_List_1_2025-10-06.csv")
labels <- as.numeric(format(x = seq((min(enrichment_kegg$FDR ) - 1e-8), (max(enrichment_kegg$FDR) + 1e-4), length.out = 10), scientific = TRUE, digits = 3))
min(enrichment_kegg$Fold.Enrichment); max(enrichment_kegg$Fold.Enrichment)
ggplot(data = enrichment_kegg, mapping = aes(x = `Fold.Enrichment`, y = Term, colour = FDR, size = Count)) + 
    scale_x_continuous(limits = c(0, 5)) + 
    # geom_point(alpha = 0.5) + 
    # scale_size_continuous(range = c(min(5), max(15))) + 
    geom_segment(mapping = aes(x = 0, y = Term, xend = `Fold.Enrichment`, yend = Term), size = 1) + 
    geom_point(mapping = aes(size = Count)) + 
    scale_color_gradient(low = "gold", high = "red", limits = c(min(labels), max(labels)), 
                         breaks = labels,
                         labels = as.character(labels),
                         guide = guide_colorbar(
                             title = "FDR", 
                             title.position = "top", 
                             title.hjust = 0.5, 
                             barwidth = unit(1, "lines"), 
                             barheight = unit(15, "lines"), 
                             ticks = TRUE,
                             ticks.colour = "black"
                         ),
                         oob = scales::squish) + 
    theme_bw() + 
    theme(
        axis.text.x = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(t = 15)),
        axis.title.y = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(r = 15)),
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(size = 1),
        legend.title = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(b = 15)),
        legend.text = element_text(size = 15, family = "serif", face = "bold", colour = "black"),
        plot.title = element_text(size = 30, family = "serif", face = "bold", colour = "black"),
        panel.border = element_rect(linewidth = 1)
    ) + 
    labs(
        title = "KEGG Pathway Gene Enrichment Analysis of DEGs",
        x = "Fold Enrichment",
        y = "Significantly Enriched KEGG Pathways"
    )





# Lollipop Graph for Shing GO Output
enrichment_kegg <- read.csv(file = "C:\\Users\\HP\\Documents\\R_Go\\Raw Combined + Filtered\\MT-Genes Filtered Combined\\Latest_1.01_Low_Resolution_C15\\_Downstream_Analysis\\Gene Ontology Report\\GO_Threshold_0.01\\Cellular Component/enrichment.csv")

# Lollipop Chart

labels <- as.numeric(format(x = seq((min(enrichment_kegg$Enrichment.FDR ) - 1e-8), (max(enrichment_kegg$Enrichment.FDR) + 1e-4), length.out = 10), scientific = TRUE, digits = 3))
min(enrichment_kegg$Fold.Enrichment); max(enrichment_kegg$Fold.Enrichment)
ggplot(data = enrichment_kegg, mapping = aes(x = `Fold.Enrichment`, y = Pathway, colour = Enrichment.FDR, size = nGenes)) + 
    scale_x_continuous(limits = c(0, 10)) + 
    # geom_point(alpha = 0.5) + 
    # scale_size_continuous(range = c(min(5), max(15))) + 
    geom_segment(mapping = aes(x = 0, y = Pathway, xend = `Fold.Enrichment`, yend = Pathway), size = 1) + 
    geom_point(mapping = aes(size = nGenes)) + 
    scale_color_gradient(low = "gold", high = "red", limits = c(min(labels), max(labels)), 
                         breaks = labels,
                         labels = as.character(labels),
                         guide = guide_colorbar(
                             # title = "FDR", 
                             title.position = "top", 
                             title.hjust = 0.5, 
                             barwidth = unit(1, "lines"), 
                             barheight = unit(20, "lines"), 
                             ticks = TRUE,
                             ticks.colour = "black"
                         ),
                         oob = scales::squish) + 
    # scale_fill_gradient(limits = c(min(enrichment_kegg$Enrichment.FDR), max(enrichment_kegg$Enrichment.FDR)), low = "gold", high = "red") + 
    # scale_color_continuous(limits = c(min(enrichment_kegg$Enrichment.FDR), max(enrichment_kegg$Enrichment.FDR)))
    theme_bw() + 
    theme(
        axis.text.x = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(t = 15)),
        axis.title.y = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(r = 15)),
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(size = 1),
        legend.title = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(b = 15)),
        legend.text = element_text(size = 15, family = "serif", face = "bold", colour = "black"),
        plot.title = element_text(size = 30, family = "serif", face = "bold", colour = "black"),
        panel.border = element_rect(linewidth = 1)
    ) + 
    labs(
        title = "Gene Ontology Analysis of DEGs",
        x = "Fold Enrichment",
        y = "Significantly Enriched Cellular Components"
    )
# Good









min(enrichment_kegg$Fold.Enrichment); max(enrichment_kegg$Fold.Enrichment)
# For Bubble Chart
ggplot(data = enrichment_kegg, mapping = aes(x = Fold.Enrichment, y = Pathway, colour = Enrichment.FDR, size = nGenes)) + 
    scale_x_continuous(limits = c(0, 15)) + 
    geom_point(alpha = 0.5) + 
    scale_size_continuous(range = c(min(5), max(15))) + 
    # geom_segment(mapping = aes(x = 0, y = Term, xend = `Fold Enrichment`, yend = Term), size = 1) + 
    # geom_point(mapping = aes(size = Count)) + 
    scale_color_gradient(low = "gold", high = "red", limits = c(min(labels), max(labels)), 
                         breaks = labels,
                         labels = as.character(labels),
                         guide = guide_colorbar(
                             title = "FDR", 
                             title.position = "top", 
                             title.hjust = 0.5, 
                             barwidth = unit(1, "lines"), 
                             barheight = unit(15, "lines"), 
                             ticks = TRUE,
                             ticks.colour = "black"
                         ),
                         oob = scales::squish) + 
    theme_bw() + 
    theme(
        axis.text.x = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(t = 15)),
        axis.title.y = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(r = 15)),
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 1),
        legend.title = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(b = 15)),
        legend.text = element_text(size = 15, family = "serif", face = "bold", colour = "black"),
        plot.title = element_text(size = 30, family = "serif", face = "bold", colour = "black"),
        panel.border = element_rect(linewidth = 1)
    ) + 
    labs(
        title = "KEGG Pathway Gene Enrichment Analysis",
        x = "Fold Enrichment",
        y = "Significantly Enriched KEGG Pathways"
    )
