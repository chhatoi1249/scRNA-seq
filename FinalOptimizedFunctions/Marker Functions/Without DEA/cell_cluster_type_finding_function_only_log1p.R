find_cell_cluster_types <- function(mat, cell_markers, clusters){
    temp_mat_log1p <- log1p(mat)
    
    # clusters: Holds all cluster information.
    # cell_markers: Holds all cell types identified!
    # mat_log1p: Log base exp(1) normalized sparse matrix
    
    type_markers <- vector()
    
    for (i in levels(clusters)) {
        mat_log1p <- as.matrix(temp_mat_log1p[,clusters == i])
        
        l <- data.frame(X = colnames(mat_log1p))
        for (i in 1:length(cell_markers)){
            common_genes <- intersect(cell_markers[[i]], rownames(mat_log1p))
            if(length(common_genes) == 0){
                temp <- rep(NA, ncol(mat_log1p))
            }else{
                temp <- colMeans(mat_log1p[common_genes, ,drop = F])
            }
            l <- cbind(l, as.data.frame(temp))
        }
        l <- l[-1]
        names(cell_markers) -> names(l)
        
        colMeans_DataFrame <- colMeans(l)
        
        type_markers <- append(type_markers, names(which.max(colMeans_DataFrame)))
        
        if (is.null(names(which.max(colMeans_DataFrame)))) {
            type_markers <- append(type_markers, NA)
        }
    }
    
    return(type_markers)
}
