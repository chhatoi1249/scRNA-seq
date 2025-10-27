filter_library_size <- function(n = 1){
    message("If you're not acquainted with a specific field parameter in this process, simply leave it blank.\n\nMake sure you entered the correct path for your file(s).\n\nThese file(s) have a .csv extension.\n\nMake sure your CSV file's first row contains all the barcodes and the first column contains all the gene IDs.\n\nYou'll get a list named as 'lSize_Matrix' after this process was successfully executed, which holds the filtered matrix (matrices) by library size.")
    # n <- readline(prompt = "Total No. of files: ")
    n <- n
    if(n != ""){
        n <- as.integer(n)
    }
    lSize_Matrix <- list()
    if(is.integer(n) & n > 0){
        path_vec <- vector()
        vec_k <- vector()
        for(i in 1:n){
            path <- readline(prompt = "Please enter the file path for your count data CSV file.\nPath:")
            k <- readline(prompt = "Threshold library size for this data: ")
            if(k != "" & path != ""){
                k <- as.integer(k)
                vec_k <- append(vec_k, k)
                path_vec <- append(path_vec, path)
            }
        }
        rm(path, k, n, i)
        if(length(vec_k) != 0){
            while(length(vec_k) != 0){
                t <- read.csv(path_vec[1], row.names = 1) |> as.matrix()
                lSize_Matrix[[length(lSize_Matrix)+1]] <- t[,(colSums(t) >= vec_k[1])] |> as.data.frame()
                vec_k <- vec_k[-1];path_vec <- path_vec[-1]
            }
            rm(path_vec, vec_k, t)
        }
    }
    
    assign("lSize_Matrix",lSize_Matrix, envir = as.environment(1))
    message("This process has been completed successfully.")
}
