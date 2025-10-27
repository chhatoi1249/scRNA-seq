# CapEff Plot
x <- as.data.frame(capEff)
ggplot(data = x, aes(x = capEff)) + geom_histogram(bins = 2000, colour = "darkred", fill = "red", alpha = .5) + 
    theme_bw() + 
    geom_density(colour = "black", fill = "gold", alpha = .6)

ggplot(data = x, aes(x = capEff)) + geom_density(colour = "black", fill = "red", alpha = .5) + 
    theme_bw()

# Volcano Plot


# Prepare the dataset


# Volcano Plot
# SwarnAdjLRT_Dataset <- as.data.frame(SwarnAdjLRT_4030_Clusters_5_DEseq.norm$SwarnAdjLRT)

SwarnAdjLRT_4030_Clusters_15.DEseq.norm <- SwarnAdjLRT_500_Iter

SwarnAdjLRT_Dataset <- as.data.frame(SwarnAdjLRT_4030_Clusters_15.DEseq.norm$SwarnAdjLRT)

# SwarnAdjLRT_Dataset <- SwarnAdjLRT_Dataset |> filter(DE.Adj.Pvalue > 0)

# SwarnAdjLRT_Dataset$DE.Adj.Pvalue <- -log10(SwarnAdjLRT_Dataset$DE.Adj.Pvalue)


max_value <- max(SwarnAdjLRT_Dataset$Adj_Log2_Fold_Change[!is.infinite(SwarnAdjLRT_Dataset$Adj_Log2_Fold_Change)])

min_value <- min(SwarnAdjLRT_Dataset$Adj_Log2_Fold_Change[!is.infinite(SwarnAdjLRT_Dataset$Adj_Log2_Fold_Change)])

require(dplyr); require(EnhancedVolcano); require(ggplot2)

SwarnAdjLRT_Dataset <- SwarnAdjLRT_Dataset |> mutate(Adj_Log2_Fold_Change = if_else(is.infinite(Adj_Log2_Fold_Change), max_value * 1.5, Adj_Log2_Fold_Change))

SwarnAdjLRT_Dataset <- SwarnAdjLRT_Dataset |> mutate(Adj_Log2_Fold_Change = if_else(Adj_Log2_Fold_Change == 0, min_value * 1.5, Adj_Log2_Fold_Change))



# SwarnAdjLRT_Dataset <- SwarnAdjLRT_Dataset |> mutate(DE.Adj.Pvalue = if_else(DE.Adj.Pvalue == 0, 1E50, DE.Adj.Pvalue))


# max(SwarnAdjLRT_Dataset$DE.Adj.Pvalue[!is.infinite(SwarnAdjLRT_Dataset$DE.Adj.Pvalue)])
# min(SwarnAdjLRT_Dataset$DE.Adj.Pvalue[!is.infinite(SwarnAdjLRT_Dataset$DE.Adj.Pvalue)])
# -log10(min(SwarnAdjLRT_Dataset$DE.Adj.Pvalue[!is.infinite(SwarnAdjLRT_Dataset$DE.Adj.Pvalue)]))

dataset <- SwarnAdjLRT_Dataset |> select(Adj_Log2_Fold_Change, DE.Adj.Pvalue)
# save(dataset, file = "VolcanoPlot.csv")
# dataset <- dataset[dataset$DE.Adj.Pvalue > 0, ] # Never skip because the are highly significant for the study. 

# dataset <- cbind(dataset, log2FC = log2FC)
# EnhancedVolcano(toptable = dataset, lab = rownames(dataset), x = "log2FC", y = "DE.Adj.Pvalue", pCutoff = 1e-5, ylim = c(-6, 350), xlim = c(-7, 7))

dataset <- dataset |> mutate(DE.Adj.Pvalue = if_else(dataset$DE.Adj.Pvalue == 0, min(dataset$DE.Adj.Pvalue[dataset$DE.Adj.Pvalue!=0]) * 10^-3, dataset$DE.Adj.Pvalue))

EnhancedVolcano(toptable = dataset, lab = rownames(dataset), x = "Adj_Log2_Fold_Change", y = "DE.Adj.Pvalue", pCutoff = 1e-5, ylim = c(-5, 350), xlim = c(-5, 5))

# dataset$DE.Adj.Pvalue[dataset$DE.Adj.Pvalue == 10^-320] -> m

new_data_i <- SwarnAdjLRT_Dataset |> filter(Adj_Log2_Fold_Change > 1 & DE.Adj.Pvalue < 1e-5)
dim(new_data_i)
i_names <- rownames(new_data_i)

new_data_c <- SwarnAdjLRT_Dataset |> filter(Adj_Log2_Fold_Change < -1 & DE.Adj.Pvalue < 1e-5)
dim(new_data_c)
c_names <- rownames(new_data_c)



DEGenes <- c(i_names, c_names)

DEGData <- ExtDEGAdjNormData(CountData = combined_filtered_matrix, DEGenes = DEGenes, norm.method = "DEseq.norm", group = v, RNAspike.use = F)


control <- DEGData$control
infected <- DEGData$infected

dim(control); dim(infected)

length(c_names); length(i_names)


data <- cbind(DEGData$control, DEGData$infected)
