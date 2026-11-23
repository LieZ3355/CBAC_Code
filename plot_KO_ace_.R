plot_KO_ace_ <- function (X, gKO, q = 0.99, annotate = TRUE, nCategories = 20, 
          fdrThreshold = 0.05) 
{
  gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$distance > 
                                                 1e-10 & X$diffRegulation$p.adj < 0.05]))
  if (length(gList) > 0) {
    sCluster <- as.matrix(X$tensorNetworks$WT[gList, gList])
    koInfo <- sCluster[gKO, ]
    gList <- gList[!grepl("^mt-|^Rpl|^Rps", gList, ignore.case = TRUE)]
    sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
    sCluster[gKO, ] <- koInfo
    diag(sCluster) <- 0
    sCluster <- reshape2::melt(as.matrix(sCluster))
    colnames(sCluster) <- c("from", "to", "W")
    sCluster <- sCluster[sCluster$W != 0, ]
    netPlot <- igraph::graph_from_data_frame(sCluster, directed = TRUE)
    dPlot <- igraph::centr_degree(netPlot)$res
    W <- rep(1, nrow(sCluster))
    sG <- (names(igraph::V(netPlot))[dPlot > 1])[-1]
    W[sCluster$from %in% sG] <- 0.2
    W[sCluster$to %in% sG] <- 0.2
    W[sCluster$from %in% gKO] <- 1
    W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
    set.seed(1)
    layPlot <- igraph::layout_with_fr(netPlot, weights = W)
    dPlot <- (dPlot/max(dPlot)) * 20
    if (isTRUE(annotate)) {
      enrichFunction <- function(X, fdrThreshold = fdrThreshold) {
        En_data <- enrichr(X, c("KEGG_2019_Human", "GO_Biological_Process_2018", 
                          "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", 
                          "BioPlanet_2019", "WikiPathways_2019_Human", 
                          "Reactome_2016"))
        En_data <- do.call(rbind.data.frame, En_data)
        En_data <- En_data[En_data$Adjusted.P.value < fdrThreshold, ]
        En_data <- En_data[order(En_data$Adjusted.P.value), ]
        En_data$Term <- unlist(lapply(strsplit(En_data$Term, ""), 
                                function(X) {
                                  X[1] <- toupper(X[1])
                                  X <- paste0(X, collapse = "")
                                  X <- gsub("\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+", 
                                            "", X)
                                  X <- gsub("'s", "", X)
                                  X <- unlist(strsplit(X, ","))[1]
                                  X <- gsub("[[:blank:]]$", "", X)
                                  return(X)
                                }))
        selectedSet <- rep(FALSE, nrow(En_data))
        for (i in seq_len(nrow(En_data))) {
          if (i == 1) {
            selectedSet[i] <- TRUE
          }
          else {
            A <- unique(unlist(strsplit(En_data[which(selectedSet[seq_len(i)]), 
                                          "Genes"], ";")))
            B <- unlist(strsplit(En_data[i, "Genes"], ";"))
            selectedSet[i] <- !all(B %in% A)
          }
        }
        gSets <- table(toupper(En_data$Term))
        gSets <- names(gSets[gSets > 1])
        for (i in gSets) {
          selectedSet[which(toupper(En_data$Term) %in% i)[-1]] <- FALSE
        }
        En_data <- En_data[selectedSet, ]
        if (nrow(En_data) > nCategories) {
          En_data <- En_data[seq_len(nCategories), ]
        }
        return(En_data)
      }
      En_data <- enrichFunction(gList, fdrThreshold)
      if (isTRUE(nrow(En_data) > 0)) {
        tPlot <- strsplit(En_data$Genes, ";")
        pPlot <- matrix(0, nrow = length(igraph::V(netPlot)), 
                        ncol = nrow(En_data))
        rownames(pPlot) <- toupper(names(igraph::V(netPlot)))
        for (i in seq_along(tPlot)) {
          pPlot[unlist(tPlot[i]), i] <- 1
        }
        pPlot <- lapply(seq_len(nrow(pPlot)), function(X) {
          as.vector(pPlot[X, ])
        })
        names(pPlot) <- names(igraph::V(netPlot))
        tPlot <- unique(unlist(tPlot))
        eGenes <- toupper(names(igraph::V(netPlot))) %in% 
          tPlot
        vColor <- rgb(195/255, 199/255, 198/255, 0.3)
        if (nrow(En_data) == 1) {
          pieColors <- list(hcl.colors(5, palette = "Zissou 1", 
                                       alpha = 0.7)[5])
        }
        else {
          pieColors <- list(hcl.colors(nrow(En_data), palette = "Zissou 1", 
                                       alpha = 0.7))
        }
        par(mar = c(4, 0, 0, 0), xpd = TRUE)
        suppressWarnings(plot(netPlot, layout = layPlot, 
                              edge.arrow.size = 0.2, vertex.label.color = "black", 
                              vertex.shape = ifelse(eGenes, "pie", "circle"), 
                              vertex.pie = pPlot, vertex.size = 10 + dPlot, 
                              vertex.pie.color = pieColors, vertex.label.family = "Times", 
                              vertex.label.font = ifelse(eGenes, 2, 1), 
                              edge.color = ifelse(igraph::E(netPlot)$W > 0, "red", 
                                                  "blue"), edge.curved = ifelse(W == 0.2, 
                                                                                0, 0.1), vertex.color = vColor, vertex.frame.color = NA))
        sigLevel <- formatC(En_data$Adjusted.P.value, digits = 2, 
                            format = "g", width = 0, drop0trailing = TRUE)
        gSetNames <- lengths(strsplit(En_data$Genes, ";"))
        gSetNames <- paste0("(", gSetNames, ") ", En_data$Term, 
                            " FDR = ", sigLevel)
        legend(x = -1.05, y = -1.05, legend = gSetNames, 
               bty = "n", ncol = 2, cex = 1, col = unlist(pieColors), 
               pch = 16)
      }
      else {
        par(mar = c(0, 0, 0, 0))
        plot(netPlot, layout = layPlot, edge.arrow.size = 0.2, 
             vertex.label.color = "black", vertex.size = 10 + 
               dPlot, vertex.label.family = "Times", edge.color = ifelse(igraph::E(netPlot)$W > 
                                                                           0, "red", "blue"), edge.curved = ifelse(W == 
                                                                                                                     0.2, 0, 0.1), vertex.color = rgb(0, 188/255, 
                                                                                                                                                      1, 0.3), vertex.frame.color = NA)
      }
    }
    else {
      par(mar = c(0, 0, 0, 0))
      plot(netPlot, layout = layPlot, edge.arrow.size = 0.2, 
           vertex.label.color = "black", vertex.size = 10 + 
             dPlot, vertex.label.family = "Times", edge.color = ifelse(igraph::E(netPlot)$W > 
                                                                         0, "red", "blue"), edge.curved = ifelse(W == 
                                                                                                                   0.2, 0, 0.1), vertex.color = rgb(0, 188/255, 
                                                                                                                                                    1, 0.3), vertex.frame.color = NA)
    }
  }
  else {
    netPlot <- matrix(0, 1, 1)
    rownames(netPlot) <- colnames(netPlot) <- gKO
    netPlot <- igraph::graph_from_adjacency_matrix(netPlot)
    par(mar = c(0, 0, 0, 0))
    plot(netPlot, edge.arrow.size = 0.2, vertex.size = 50, 
         vertex.label.color = "black", vertex.label.family = "Times", 
         vertex.color = rgb(0, 188/255, 1, 0.3), vertex.frame.color = NA)
  }
}
