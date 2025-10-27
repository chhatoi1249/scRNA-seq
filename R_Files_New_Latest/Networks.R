fun <- function(linklist) {
    x <- linklist
    x$regulatoryGene <- as.character(x$regulatoryGene)
    x$targetGene <- as.character(x$targetGene)
    vec <- c()
    # not_vec <- c()
    for (i in 1:dim(x)[1]) {
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

###################################################################################################################################
###################################################################################################################################
# Control samples
control_linklist <- GENIE3::getLinkList(weightMatrix = GENIE3_1.09_Control_DEGenes, threshold = .10)

control_linklist <- fun(linklist = control_linklist)

library(igraph)
control_graph <- igraph::graph_from_data_frame(d = control_linklist, directed = TRUE)

##############################################################################
V(control_graph)$degree <- degree(control_graph, mode = "out")
# Calculate out-degree (number of targets regulated)
control_out_degree <- degree(control_graph, mode = "out")

# Or calculate weighted out-degree (sum of outgoing edge weights)
control_weighted_out_degree <- strength(control_graph, mode = "out", weights = E(control_graph)$weight)
##############################################################################


V(control_graph)$color = "red";
E(control_graph)$color = "black";

# plot(graph,vertex.size = V(graph)$degree*.35,edge.arrow.size = .4,layout = layout.sphere,vertex.size = c(20:40),vertex.shape = "circle",vertex.frame.color = "black",vertex.size2 = 40,vertex.label.color = "black",vertex.label.font = 4,vertex.label.cex = 1.2,vertex.label.dist = .5,vertex.label.degree = 20,edge.curved = 1)

# layout = layout.davidson.harel

plot(control_graph,vertex.size = control_weighted_out_degree*20,edge.arrow.size = .4,layout = layout.sphere,vertex.size = c(20:40),vertex.shape = "circle",vertex.frame.color = "black",vertex.size2 = 40,vertex.label.color = "black",vertex.label.font = 4,vertex.label.cex = 1.5,vertex.label.dist = .5,vertex.label.degree = 20,edge.curved = 1, fontsize = 2)

# Find top hub genes by weighted out-degree
control_top_hubs <- sort(control_weighted_out_degree, decreasing = TRUE)
# control_top_hubs <- top_hubs[control_top_hubs > 0]
length(control_top_hubs)
control_top_hubs[control_top_hubs > 0] |> length()
###################################################################################################################################
###################################################################################################################################

# Suppose gene A regulates 3 genes with weights 0.8, 0.5, and 0.3.
# Out-degree of gene A = 3 (three outgoing edges)
# Weighted out-degree of gene A = 0.8 + 0.5 + 0.3 = 1.6
# 
# 
# Gene B regulates 5 genes but with weights 0.1 each.+
# Out-degree of gene B = 5
# Weighted out-degree of gene B = 0.1 * 5 = 0.5
# Even though gene B regulates more genes, gene A has a higher weighted out-degree, indicating stronger overall regulatory influence.






###################################################################################################################################
###################################################################################################################################
# Infected samples
GENIE3::GENIE3(exprMatrix = )
infected_linklist <- GENIE3::getLinkList(weightMatrix = GENIE3_1.09_Infected_DEGenes, threshold = .10)

infected_linklist <- fun(linklist = infected_linklist)

library(igraph)
infected_graph <- igraph::graph_from_data_frame(d = infected_linklist, directed = TRUE)


##############################################################################
V(infected_graph)$degree <- degree(infected_graph, mode = "out")
# Calculate out-degree (number of targets regulated)
infected_out_degree <- degree(infected_graph, mode = "out")

# Or calculate weighted out-degree (sum of outgoing edge weights)
infected_weighted_out_degree <- strength(infected_graph, mode = "out", weights = E(infected_graph)$weight)
##############################################################################


V(infected_graph)$color = "red";
E(infected_graph)$color = "black";

# plot(graph,vertex.size = V(graph)$degree*.35,edge.arrow.size = .4,layout = layout.sphere,vertex.size = c(20:40),vertex.shape = "circle",vertex.frame.color = "black",vertex.size2 = 40,vertex.label.color = "black",vertex.label.font = 4,vertex.label.cex = 1.2,vertex.label.dist = .5,vertex.label.degree = 20,edge.curved = 1)

plot(infected_graph,vertex.size = infected_weighted_out_degree*20,edge.arrow.size = .4,layout = layout.davidson.harel,vertex.size = c(20:40),vertex.shape = "circle",vertex.frame.color = "black",vertex.size2 = 40,vertex.label.color = "black",vertex.label.font = 4,vertex.label.cex = 2,vertex.label.dist = .5,vertex.label.degree = 20,edge.curved = 1)


# Find top hub genes by weighted out-degree
infected_top_hubs <- sort(infected_weighted_out_degree, decreasing = TRUE)
length(infected_top_hubs)
length(infected_top_hubs[infected_top_hubs > 0])
###################################################################################################################################
###################################################################################################################################



















###################################################################################################################################
###################################################################################################################################

# Compare control Vs infected hub genes
control_hub_genes <- names(control_top_hubs[control_top_hubs > 0]) #control_hub_genes[1:10] # 
infected_hub_genes <- names(infected_top_hubs[infected_top_hubs > 0])
# infected_hub_genes[1:10] # 
m <- intersect(control_hub_genes, infected_hub_genes)

write.csv(x = control_hub_genes[!(control_hub_genes %in% m)], row.names = FALSE, file = "Uninfected_0.10_Hub_Genes.csv")

write.csv(x = infected_hub_genes[!(infected_hub_genes %in% m)], row.names = FALSE, file = "H5N1_Infected_0.10_Hub_Genes.csv")

write.csv(m, file = "HK_0.10_Genes.csv")


intersect(control_hub_genes, infected_hub_genes) |> length()
`Uninfected Hub Genes` <- control_hub_genes
`H5N1 Infected Hub Genes` <- infected_hub_genes
list <- list(`Uninfected Hub Genes` = `Uninfected Hub Genes`, `H5N1 Infected Hub Genes` = `H5N1 Infected Hub Genes`)
par(mar = c(6, 8, 4, 2))
UpSetR::upset(data = UpSetR::fromList(list), nsets = 2, nintersects = NA, keep.order = TRUE, matrix.color = "#AA12C8", main.bar.color = "#30BF91", mainbar.y.label = "Common Hub Genes", mainbar.y.max = 400, sets.bar.color = "red", sets.x.label = "Number of Hub Genes \nIdentified per Case", point.size = 30, line.size = 6, mb.ratio = c(.6, .4), order.by = "freq", number.angles = 0, shade.color = "#E2E0BD", matrix.dot.alpha = .5, set_size.show = TRUE, set_size.numbers_size = 12, set_size.angles = 0, text.scale = c(6, 5, 6, 5, 8, 6), set_size.scale_max = 550)


col <- c(Set1 = "#1f78b4", Set2 = "#33a02c", Set3 = "#ff7f00")

###################################################################################################################################
###################################################################################################################################
































###################################################################################################################################
###################################################################################################################################
# All combined analysis (Optional)
# Combined samples
combined_linklist <- GENIE3::getLinkList(weightMatrix = genie3_results_only_degenes_data_normalized, threshold = .1)

combined_linklist <- fun(linklist = combined_linklist)

library(igraph)
combined_graph <- igraph::graph_from_data_frame(d = combined_linklist, directed = TRUE)


##############################################################################
V(combined_graph)$degree <- degree(combined_graph, mode = "out")
# Calculate out-degree (number of targets regulated)
combined_out_degree <- degree(combined_graph, mode = "out")

# Or calculate weighted out-degree (sum of outgoing edge weights)
combined_weighted_out_degree <- strength(combined_graph, mode = "out", weights = E(combined_graph)$weight)
##############################################################################


# Cut this paragraph
# unique(c(combined_linklist$regulatoryGene, combined_linklist$targetGene)) -> combined_unique_all_vertices
# combined_graph <- induced_subgraph(graph = combined_graph, vids = unique(combined_unique_all_vertices %in% combined_linklist$regulatoryGene))



V(combined_graph)$color = "red";
E(combined_graph)$color = "black";

# plot(graph,vertex.size = V(graph)$degree*.35,edge.arrow.size = .4,layout = layout.sphere,vertex.size = c(20:40),vertex.shape = "circle",vertex.frame.color = "black",vertex.size2 = 40,vertex.label.color = "black",vertex.label.font = 4,vertex.label.cex = 1.2,vertex.label.dist = .5,vertex.label.degree = 20,edge.curved = 1)

plot(combined_graph,vertex.size = combined_weighted_out_degree*30,edge.arrow.size = .4,layout = layout.sphere,vertex.size = c(20:40),vertex.shape = "circle",vertex.frame.color = "black",vertex.size2 = 40,vertex.label.color = "black",vertex.label.font = 4,vertex.label.cex = 1.2,vertex.label.dist = .5,vertex.label.degree = 20,edge.curved = 1)


# Find top hub genes by weighted out-degree
combined_top_hubs <- sort(combined_weighted_out_degree, decreasing = TRUE)
length(combined_top_hubs)
length(combined_top_hubs[combined_top_hubs > 0])


# Optional Upset plot
# Compare control Vs infected hub genes
control_hub_genes <- names(control_top_hubs[control_top_hubs > 0]) #control_hub_genes[1:10] # 
infected_hub_genes <- names(infected_top_hubs[infected_top_hubs > 0])
# infected_hub_genes[1:10] #
combined_hub_genes <- names(combined_top_hubs[combined_top_hubs > 0])

intersect(control_hub_genes, infected_hub_genes) |> length()
list <- list(control_hub_genes = control_hub_genes, infected_hub_genes = infected_hub_genes, combined_hub_genes = combined_hub_genes)

UpSetR::upset(data = UpSetR::fromList(list), nsets = 3, nintersects = NA, keep.order = TRUE, matrix.color = "#AA12C8", main.bar.color = "#FD5E7B", mainbar.y.label = "No. of Common Hub Genes", mainbar.y.max = 400, sets.bar.color = "red", sets.x.label = "Number of Genes in Individual Sets", point.size = 30, line.size = 6, mb.ratio = c(.7, .3), order.by = "freq", number.angles = 0, shade.color = "#E2E0BD", matrix.dot.alpha = .5, set_size.show = TRUE, set_size.numbers_size = 12, set_size.angles = 0, text.scale = c(5, 4, 5, 4, 6, 7), set_size.scale_max = 500)
###################################################################################################################################
###################################################################################################################################