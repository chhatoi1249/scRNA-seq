Filter_MT_Content <- function(counts_data, MT_Threshold) {
    names <- rownames(counts_data)
    vec <- vector()
    for (i in 1:length(names)) {
        if (startsWith(names[i], "MT")) {
            vec <- append(vec, names[i])
        }
    }
    
    org_cs <- colSums(counts_data)
    mt_cs <- colSums(counts_data[vec,])
    out <- mt_cs/org_cs
    
    counts_data <- counts_data[,out <= MT_Threshold]
    
    return(counts_data)
}
