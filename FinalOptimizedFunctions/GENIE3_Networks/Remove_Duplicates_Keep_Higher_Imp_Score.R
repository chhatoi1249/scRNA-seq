fun <- function(linklist) {
    x <- linklist
    x$regulatoryGene <- as.character(x$regulatoryGene)
    x$targetGene <- as.character(x$targetGene)
    vec <- c()
    # not_vec <- c()
    for (i in seq_len(length.out = dim(x)[1])) {
        if (x$regulatoryGene[i] %in% x$targetGene) {
            genes <- which(x$targetGene %in% x$regulatoryGene[i] == TRUE)
            for (j in genes) {
                if (x$targetGene[i] == x$regulatoryGene[j]) {
                    if (x$weight[i] >= x$weight[j]) {
                        vec <- append(x = vec, values = j)
                        # not_vec <- append(x = not_vec, values = i)
                    }else {
                        vec <- append(x = vec, values = i)
                        # not_vec <- append(x = not_vec, values = j)
                    }
                }
            }
        }
    }
    if (is.null(vec)) return(x) else return(x[-vec, ])
}