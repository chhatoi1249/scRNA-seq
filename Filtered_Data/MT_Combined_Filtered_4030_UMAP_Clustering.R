
########################## Filter the dataset ########################

filter_gene_IDs()
write.csv(Result_Matrix[[1]], file = "Combined_Gene_IDs_Filtered_Matrix.csv")

filter_library_size()
write.csv(lSize_Matrix[[1]], file = "Combined_Final_Filtered_Matrix.csv")

############################# Load required packages ####################

require(Matrix); require(Seurat); require(MultiK); require(dplyr); require(ggplot2); require(ggrepel)

############################# Import data ###############################

combined_filtered_matrix <- read.csv("Final_Combined_Filtered_Matrix.csv", row.names = 1) |> as.matrix() |> Matrix(sparse = T)

Filter_MT_Content(counts_data = combined_filtered_matrix, MT_Threshold = .05) -> combined_filtered_matrix

############################# Create seurat object ############################

seu_combined_object <- CreateSeuratObject(counts = combined_filtered_matrix, project = "MT Content Filtered Combined sc-RNA-seq Analysis")

############################# Preprocessing ##########################

seu_combined_object <- NormalizeData(seu_combined_object) |> FindVariableFeatures() |> ScaleData() |> RunPCA()

############################# Run MultiK ############################

filtered_combined_multik_results <- MultiK(seu = seu_combined_object, resolution = seq(0.05, 2, .01), nPC = 30, reps = 100, pSample = .8)

############################ MultiK Graph ################################

DiagMultiKPlot(filtered_combined_multik_results)

############################# Optimal K values ###############################

# Optimal K: 15 11  5
# Corresponding resolutions: 1.11 .75 .18

############################### UMAP and find neighbours ########################

seu_combined_object <- FindNeighbors(seu_combined_object, dims = 1:50)
seu_combined_object <- FindClusters(seu_combined_object, resolution = 1.11)

seu_combined_object <- RunUMAP(seu_combined_object, dims = 1:50)

############################### Create cell markers list ################

# Constant

markers_list <- create_marker_list(excel_sheet_path = "../R_Go/Supplementary Table 1-3.xlsx")

################################ Find the cell types ##############################

# Constant

cell_types <- find_cell_types(mat_path = "Combined_Final_Filtered_Matrix.csv", cell_markers = markers_list)

################################ Now save the annotated matrix file #############################

# Constant

d <- combined_filtered_matrix |> as.data.frame()
sum(names(cell_types) == colnames(d))
names(d) <- cell_types
write.csv(d, file = "Cell_Annotated_Filtered_Combined_Matrix.csv")
rm(d)

################################ Find the markers cluster wise ################

cluster_markers <- FindAllMarkers(seu_combined_object, only.pos = T, min.pct = .25, logfc.threshold = .25)

top_cluster_markers <- cluster_markers |> group_by(cluster) |> top_n(11627, wt = avg_log2FC)

################################ Now find the cluster types ##############################

cluster_types <- find_cluster_types(cell_annotated_csv_file_path = "Cell_Annotated_Filtered_Combined_Matrix.csv", top_markers = top_cluster_markers)

# cluster_types_new <- find_cell_cluster_types(mat = combined_filtered_matrix, cell_markers = markers_list, clusters = seu_combined_object$seurat_clusters)
cluster_types_new <- find_cell_cluster_types(mat = combined_filtered_matrix, cell_markers = markers_list, clusters = seu_combined_object$seurat_clusters)

# cluster_types_new <- paste(c("Cluster_1:", "Cluster_2:", "Cluster_3:", "Cluster_4:", "Cluster_5:", "Cluster_6:", "Cluster_7:", "Cluster_8:", "Cluster_9:", "Cluster_10:", "Cluster_11:", "Cluster_12:", "Cluster_13:", "Cluster_14:", "Cluster_15:"), cluster_types_new)
# cluster_types_new <- paste(c("Cluster_1:", "Cluster_2:", "Cluster_3:", "Cluster_4:", "Cluster_5:", "Cluster_6:", "Cluster_7:", "Cluster_8:", "Cluster_9:", "Cluster_10:", "Cluster_11:"), cluster_types_new)
# cluster_types_new <- paste(c("Cluster_1:", "Cluster_2:", "Cluster_3:", "Cluster_4:", "Cluster_5:"), cluster_types_new)

names(cluster_types_new) <- levels(seu_combined_object$seurat_clusters)


# find_cell_cluster_types(mat = counts_data, cell_markers = cell_markers, clusters = clusters) -> cluster_types

################################## Add the new factor ########################

seu_combined_object$cluster_types <- factor(seu_combined_object$seurat_clusters, levels = names(cluster_types_new), labels = cluster_types_new)

################################## Final UMAP ########################## 30 * 13 ######

g <- DimPlot(seu_combined_object, group.by = "cluster_types", label = F, label.box = F, pt.size = 2, repel = F, reduction = "umap", label.size = 5, alpha = 1)

umap_coordinates <- Embeddings(seu_combined_object, "umap")

data_frame <- data.frame(umap_coordinates, cluster = seu_combined_object$cluster_types)

centroids <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = data_frame, mean)

g <- g + guides(color = guide_legend(override.aes = list(size = 20)))

g + geom_text_repel(data = centroids, aes(x = umap_1, y = umap_2, label = cluster), fontface = "bold", size = 7, box.padding = 10, point.padding = 3, segment.color = "black", segment.size = 0.5, max.overlaps = Inf, force = 3, force_pull = 5) +
    theme_bw() + 
    labs(title = "UMAP CLUSTERS", color = "Cell Types") + 
    theme(title = element_text(size = 30, face = "bold"),
          legend.text = element_text(size = 25, face = "bold"),
          legend.position = "left",
          legend.title = element_text(size = 30, face = "bold", hjust = .5),
          axis.title = element_text(size = 25, face = "bold"),
          axis.text.x = element_text(size = 20, face = "bold", margin = margin(t = 5)),
          axis.text.y = element_text(size = 20, face = "bold", margin = margin(r = 5)),
          axis.ticks = element_line(linewidth = 1))

################################# Done ##############################

################################# Remove created objects ##################################

rm(combined_filtered_matrix, filtered_combined_multik_results, seu_combined_object)

rm(centroids, cluster_markers, data_frame, g, markers_list, top_cluster_markers, umap_coordinates, cell_types, cluster_types)

##############################################################################################