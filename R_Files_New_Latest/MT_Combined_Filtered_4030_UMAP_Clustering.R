
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

# filtered_combined_multik_results <- MultiK(seu = seu_combined_object, resolution = seq(0.05, 2, .01), nPC = 30, reps = 100, pSample = .8)

############################ MultiK Graph ################################

DiagMultiKPlot(Final_Combined_Filtered_Matrix_MultiK_Results)

############################# Optimal K values ###############################

# Optimal K: 15 11  5
# Corresponding resolutions: 1.11 0.75 0.18

############################### UMAP and find neighbours ########################

seu_combined_object <- FindNeighbors(seu_combined_object, dims = 1:50)
seu_combined_object <- FindClusters(seu_combined_object, resolution = 1.01)

seu_combined_object <- RunUMAP(seu_combined_object, dims = 1:50)


# Cell Annotation

markers_list <- create_marker_list(excel_sheet_path = "Book1.xlsx")

cluster_types <- find_cell_cluster_types(mat = combined_filtered_matrix, cell_markers = markers_list, clusters = clusters_15_deseq.norm, norm.method = "DEseq.norm")
c <- cluster_types
# cluster_types <- c

c_types <- c()
for (i in 1:length(cluster_types)) {
    c_types <- append(x = c_types, values = paste0("Cluster_", i, ": "))
}
cluster_types <- paste(c_types, cluster_types)

names(cluster_types) <- levels(seu_combined_object$seurat_clusters)


# find_cell_cluster_types(mat = counts_data, cell_markers = cell_markers, clusters = clusters) -> cluster_types

################################## Add the new factor ########################

seu_combined_object$cluster_types <- factor(seu_combined_object$seurat_clusters, levels = names(cluster_types), labels = cluster_types)

################################## Final UMAP ##########################

g <- DimPlot(seu_combined_object, group.by = "cluster_types", label = F, label.box = F, pt.size = 1.2, repel = F, reduction = "umap", label.size = 7, alpha = 1)

umap_coordinates <- Embeddings(seu_combined_object, "umap")

data_frame <- data.frame(umap_coordinates, cluster = seu_combined_object$cluster_types)

centroids <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = data_frame, mean)

g + geom_text_repel(data = centroids, aes(x = umap_1, y = umap_2, label = cluster), fontface = "bold", size = 7, box.padding = 6, point.padding = 6, segment.color = "black", segment.size = 0.5, max.overlaps = Inf, force = 10, force_pull = 3, max.iter = 1000, direction = "both", min.segment.length = 0) +
    theme_bw() + 
    labs(title = "UMAP CLUSTERS", color = "Cell Types") + 
    theme(title = element_text(size = 30, face = "bold"),
          legend.text = element_text(size = 20, face = "bold"),
          legend.position = "left",
          legend.title = element_text(size = 30, face = "bold", hjust = .5),
          axis.title = element_text(size = 25, face = "bold"),
          axis.text.x = element_text(size = 20, face = "bold", margin = margin(t = 5)),
          axis.text.y = element_text(size = 20, face = "bold", margin = margin(r = 5)),
          axis.ticks = element_line(linewidth = 1)) + 
    guides(color = guide_legend(override.aes = list(size = 12, shape = 16)))

################################# Done ##############################

seu_combined_object$seurat_clusters -> clusters_15_deseq.norm
sum(clusters_15_deseq.norm == clusters_15)
