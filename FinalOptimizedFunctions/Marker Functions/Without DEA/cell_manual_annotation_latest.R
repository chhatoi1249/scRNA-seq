find_cell_cluster_types <- function(mat, cell_markers, clusters, norm.method = c("DEseq.norm", "TMM", "log1p")){
    mat <- mat[rowSums(mat) > 0,]
    mat <- mat[,colSums(mat) > 0]
    if (is.null(dim(mat))) {
        show.custom.warning("There may be an error in the input data dimensions. This could be due to having fewer than two genes with at least one read across all or some of the cells.")
        return(invisible(NULL))
    }
    if (norm.method == "DEseq.norm") {
        GM <- function(x) exp(mean(log(x[x > 0])))
        geomMean <- apply(mat, 1, GM)
        fx <- function(x) x/geomMean
        samp <- apply(mat, 2, fx)
        fxx <- function(xx) stats::median(xx[xx != 0])
        size <- apply(samp, 2, fxx)
        fxxx <- function(xxx) xxx/size
        mat <- t(apply(mat, 1, fxxx))
        rm(geomMean, samp, size, GM, fx, fxx, fxxx)
    }else if (norm.method == "TMM") {
        m <- size_factor <- edgeR::calcNormFactors(mat)
        fxxx <- function(xxx) xxx/size_factor
        mat <- t(apply(mat, 1, fxxx))
        rm(m, fxxx)
    }else {mat <- log1p(mat)}
    type_markers <- vector()
    for (i in levels(clusters)) {
        mat_new <- as.matrix(mat[,clusters == i])
        l <- data.frame(X = colnames(mat_new))
        for (i in 1:length(cell_markers)){
            common_genes <- intersect(cell_markers[[i]], rownames(mat_new))
            if(length(common_genes) == 0){
                temp <- rep(NA, ncol(mat_new))
            }else{
                temp <- colMeans(mat_new[common_genes, ,drop = F])
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
