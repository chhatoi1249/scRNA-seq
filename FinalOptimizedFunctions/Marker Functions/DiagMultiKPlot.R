DiagMultiKPlot <- function(ks) {
    ks_vals <- ks$k
    res <- ks$consensus
    
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(ggrepel))
    suppressPackageStartupMessages(library(grid))
    suppressPackageStartupMessages(library(gridExtra))
    suppressPackageStartupMessages(library(cowplot))
    
    # Tabulate frequency of each K value (only those appearing more than once)
    tog <- as.data.frame(table(ks_vals)[table(ks_vals) > 1])
    names(tog)[1] <- "ks"
    
    # Calculate rPAC scores
    pacobj <- CalcPAC(x1 = 0.1, x2 = 0.9, xvec = tog$ks, ml = res)
    tog$rpac <- pacobj$rPAC
    tog$one_minus_rpac <- 1 - tog$rpac
    
    # Frequency bar plot
    freqPlot <- ggplot(data = tog, aes(x = ks, y = Freq)) +
        geom_bar(stat = "identity", fill = "grey80") +
        geom_text(aes(label = Freq), vjust = -0.3, size = 3.5) +
        geom_hline(yintercept = 100, linetype = "dashed", color = "black") +
        scale_x_discrete("K") +
        scale_y_continuous("Number of clustering runs") +
        theme_bw(base_size = 14) +
        theme(
            axis.text = element_text(color = "black"),
            panel.background = element_rect(color = "black"),
            strip.text = element_text(size = 12),
            strip.background = element_rect(fill = "white")
        )
    
    # rPAC plot
    rpacPlot <- ggplot(data = tog, aes(x = ks, y = rpac, group = 1)) +
        geom_point(shape = 21, color = "black", fill = "black", size = 2) +
        geom_line() +
        scale_x_discrete("K") +
        scale_y_continuous("rPAC") +
        theme_bw(base_size = 14) +
        theme(
            axis.text = element_text(color = "black"),
            panel.background = element_rect(color = "black"),
            strip.text = element_text(size = 12),
            strip.background = element_rect(fill = "white")
        )
    
    # Optionally return both plots
    # plot_grid(freqPlot, rpacPlot, ncol = 2)
    
    ###########################################################################################################################################################
    optk <- function(tog){
        tog.f <- tog[tog$Freq > 100 | tog$Freq == 100, ]
        hpts <- chull(tog.f[, c("one_minus_rpac", "Freq")])
        hpts <- c(hpts, hpts[1])
        ch.df <- tog.f[hpts, ]
        df <- ch.df[, c("ks", "one_minus_rpac", "Freq")]
        colnames(df) <- c("k", "x", "y")
        b <- c()
        end_points <- c()
        for (i in 1:(nrow(df) - 1)) {
            end_points[i] <- paste(as.character(df[i, ]$k), as.character(df[(i + 
                                                                                 1), ]$k), sep = "-")
            b[i] <- (df[(i + 1), ]$y - df[i, ]$y)/(df[(i + 1), ]$x - 
                                                       df[i, ]$x)
        }
        
        lineseg.df <- data.frame(end_points = end_points, slope = b)
        lineseg.df$p1 <- do.call("rbind", strsplit(lineseg.df$end_points, 
                                                   "-"))[, 1]
        lineseg.df$p2 <- do.call("rbind", strsplit(lineseg.df$end_points, 
                                                   "-"))[, 2]
        which.k <- as.character(ch.df[which.max(ch.df$Freq), ]$ks)
        
        
        
        
        
        
        # Start by checking if all slopes for 'which.k' connections are positive
        if (all(lineseg.df[lineseg.df$p1 == which.k | lineseg.df$p2 == which.k, ]$slope > 0)) {
            
            # If all are positive, assign optK as the current which.k
            optK <- which.k
            
        } else {
            # Otherwise, find the cluster(s) where slope < 0
            tmp <- which(lineseg.df[lineseg.df$p1 == which.k | lineseg.df$p2 == which.k, ]$slope < 0)
            tmp <- lineseg.df[lineseg.df$p1 == which.k | lineseg.df$p2 == which.k, ][tmp, ]
            
            # Extract the next node (k2) that is connected but different from which.k
            which.k2 <- as.character(c(tmp$p1, tmp$p2)[which(c(tmp$p1, tmp$p2) != which.k)])
            
            # Remove connections involving which.k
            lineseg.df.sub <- lineseg.df[lineseg.df$p1 != which.k & lineseg.df$p2 != which.k, ]
            
            if (all(lineseg.df.sub[lineseg.df.sub$p1 == which.k2 | lineseg.df.sub$p2 == which.k2, ]$slope > tmp$slope)) {
                
                optK <- c(which.k, which.k2)
                
            } else {
                # Find the third level connection where slope < 0
                tmp2 <- which(lineseg.df.sub[lineseg.df.sub$p1 == which.k2 | lineseg.df.sub$p2 == which.k2, ]$slope < 0)
                tmp2 <- lineseg.df.sub[lineseg.df.sub$p1 == which.k2 | lineseg.df.sub$p2 == which.k2, ][tmp2, ]
                
                # Extract which.k3 by excluding which.k and which.k2
                which.k3 <- as.character(c(tmp2$p1, tmp2$p2)[which(
                    c(tmp2$p1, tmp2$p2) != which.k & c(tmp2$p1, tmp2$p2) != which.k2
                )])
                
                # Remove connections involving which.k and which.k2
                lineseg.df.sub <- lineseg.df[
                    lineseg.df$p1 != which.k & lineseg.df$p2 != which.k &
                        lineseg.df$p1 != which.k2 & lineseg.df$p2 != which.k2, ]
                
                if (all(lineseg.df.sub[lineseg.df.sub$p1 == which.k3 | lineseg.df.sub$p2 == which.k3, ]$slope > tmp2$slope)) {
                    optK <- c(which.k, which.k2, which.k3)
                } else {
                    optK <- c(which.k, which.k2, which.k3)
                }
            }
        }
        return(optK)
    }
    ###########################################################################################################################################################
    optK <- optk(tog)
    cat("Optimal K: ", optK, "\n")
    
    scatPlot <- ggplot(data = tog, aes(x = one_minus_rpac, y = Freq)) +
        geom_point(shape = 21, color = "black", fill = "black", size = 1.5) +
        geom_path(color = "grey", alpha = 0.75, linetype = 2) +
        theme_bw(base_size = 14) +
        theme(
            axis.text = element_text(color = "black"),
            panel.background = element_rect(color = "black"),
            strip.text = element_text(size = 12),
            strip.background = element_rect(fill = "white")
        ) +
        scale_x_continuous("1 - rPAC") +
        scale_y_continuous("Number of clustering runs") +
        geom_hline(yintercept = 100, linetype = "dashed", color = "black") +
        geom_label_repel(
            aes(label = ks),
            segment.color = "grey50", size = 3
        ) +
        geom_path(
            data = tog[match(optk(tog), tog$ks), ]
        )
    plot_grid(freqPlot, rpacPlot, scatPlot, ncol = 3)
}
