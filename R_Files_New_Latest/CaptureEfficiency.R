SwarnSeq::CapEff(CountData = combined_filtered_matrix, CE.range = c(.01,.5), RNAspike.use = F, method = "") -> capEff

capEff <- as.data.frame(capEff)
library(ggplot2)
ggplot(data = capEff, mapping = aes(x = capEff)) + 
    geom_histogram(bins = 500,na.rm = TRUE, colour = "red", fill = "red") + 
    geom_density(alpha = 0, fill = "gold", linewidth = .7) + 
    theme_bw() + 
    theme(
        axis.text.x = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(t = 10)),
        axis.text.y = element_text(size = 15, face = "bold", family = "serif", colour = "black", margin = margin(r = 10)),
        axis.title.x = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(t = 15)),
        axis.title.y = element_text(size = 25, family = "serif", face = "bold", colour = "black", margin = margin(r = 15)),
        axis.ticks.length = unit(.5, "cm"),
        axis.ticks = element_line(linewidth = 1),
        legend.title = element_text(size = 25, family = "serif", face = "bold", colour = "black"),
        legend.text = element_text(size = 15, family = "serif", face = "bold", colour = "black"),
        plot.title = element_text(size = 30, family = "serif", face = "bold", colour = "black"),
        panel.border = element_rect(linewidth = 1)
    ) + 
    labs(
        title = "Cell's Capture Efficiency",
        y = "Frequency & Density"
    )
