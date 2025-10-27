create_marker_list <- function(excel_sheet_path){
    if(!require(readxl, quietly = T))
        install.packages("readxl")
    marker_genes_cell_types <- readxl::read_excel(excel_sheet_path) |> as.data.frame()
    
    cell_markers <- vector("list",dim(marker_genes_cell_types)[1])
    for(i in 1:dim(marker_genes_cell_types)[1]){
        temp <- unlist(strsplit(marker_genes_cell_types[i,2], split = "[,/]"))
        cell_markers[[i]][1:length(temp)] <- temp
        
    }
    names(cell_markers) <- marker_genes_cell_types[[1]]
    return(cell_markers)
}
