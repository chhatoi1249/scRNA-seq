filter_gene_IDs <- function(n = 1){
    message("If you're not acquainted with a specific field parameter in this process, simply leave it blank.\n\nMake sure you entered the correct path for your file(s).\n\nThese file(s) have a .csv extension.\n\nMake sure your CSV file's first row contains all the barcodes and the first column contains all the gene IDs.\n\nYou'll get a list named as 'Result_Matrix' after this process was successfully executed, which holds the filtered matrix (matrices) by gene IDs.")
    # n <- readline(prompt = "Total No. of files: ")
    n <- n
    if(n != ""){
        n <- as.integer(n)
    }
    Result_Matrix <- list()
    if(is.integer(n) & n > 0){
        path_vec <- vector()
        vec_k <- vector()
        for(i in 1:n){
            path <- readline(prompt = "Please enter the file path for your count data CSV file.\nPath:")
            k <- readline(prompt = "Threshold non-zero entity size for this data: ")
            if(k != "" & path != ""){
                k <- as.integer(k)
                vec_k <- append(vec_k, k)
                path_vec <- append(path_vec, path)
            }
        }
        rm(path, k, n)
        if(length(vec_k) != 0){
            while(length(vec_k) != 0){
                t <- read.csv(path_vec[1], row.names = 1) |> as.matrix()
                t_gene_IDs <- dim(t)[1]
                FUN = function(x) return(sum(x > 0) >= vec_k[1])
                c <- vector()
                for(i in 1:t_gene_IDs){
                    c <- append(c,FUN(t[i,]))
                }
                Result_Matrix[[length(Result_Matrix)+1]] <- t[c,] |> as.data.frame()
                vec_k <- vec_k[-1];path_vec <- path_vec[-1]
            }
            rm(path_vec ,vec_k, t, t_gene_IDs, FUN, c)
        }
    }
    
    assign("Result_Matrix",Result_Matrix, envir = as.environment(1))
    message("This process has been completed successfully.")
}
