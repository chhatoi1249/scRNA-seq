y <- read.csv("Biological Process//enrichment.csv") |> as.data.frame()

y2 <- round((y$nGenes/length(DEGenes))*100, 2)

names(y2) <- y$Pathway

z <- sort(y2)

# pie(x = y2, labels = paste(y$Pathway, "(", y2, "% )"), radius = 1, cex = 1.5, main = "GO Molecular Function", cex.main = 5, clockwise = F, init.angle = 110, font = 2)

pie(x = z, labels = names(z), radius = 1, cex = 1.5, main = "GO Molecular Function", cex.main = 5, clockwise = F, init.angle = 110, font = 2)


total <- length(DEGenes)

vec <- c()
for (i in 1:dim(y)[1]) {
    strsplit(y$Genes[i], split = " ") |> unlist() -> temp
    vec <- append(vec, temp)
}

significant_pathway_genes <- length(unique(vec))

temp <- c(total, significant_pathway_genes)
percentage <- round((temp / length(DEGenes)) * 100, digits = 2)
names(temp) <- c("DEGenes", "DEGenes in Gene Ontology Analysis")

pie(x = temp, labels = paste0(names(temp), "(", percentage, "%)"), radius = 1, cex = 1.5, main = "DEGenes involved in Gene Ontology (GO) Analysis", cex.main = 3, clockwise = F, init.angle = 110, font = 2, col = c("red", "darkred"))

