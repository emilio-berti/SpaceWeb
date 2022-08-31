#' @title Select MEMs
#'
#' @param webs list object obtained from local_web() function.
#' @param landscape a landscape, i.e. a matrix with xy coordinates.
#' @param output string, where to write output files. Only path and file name,
#'   with no extension.
#'
#' @return
select_mems <- function(webs, landscape, output, d2 = 0.3, style = "B") {
  # response to select MEMs
  pca <- dudi.pca(webs$properties, scannf = FALSE, nf = 2)
  if (inertia(pca)$tot.inertia[, "cum(%)"][2] < 0.7) {
    warning(
      "Cumulative variance explained by PCA = ",
      round(inertia(pca)$tot.inertia[, "cum(%)"][2]),
      "%"
    )
  }
  y <- pca$tab[, c(1, 2)]
  # candidate MEMs
  cands <- listw.candidates(
    landscape,
    d2 = d2,
    style = style,
    weights = c("flin", "fdown", "fup"),
    nb = c("del", "gab", "rel", "mst", "dnear"),
    y_fdown = seq(2, 5, by = 0.5),
    y_fup = seq(1, 5, by = 0.5)
  )
  sel <- listw.select(y, cands)
  if (is.null(sel$best.id)) {
    stop("No MEM selected. Try another style than ", style)
  }
  message("Best: ", names(sel$best.id))
  # selection summary (R2)
  r2 <- sel$candidates[order(sel$candidates$R2Adj.select, decreasing = TRUE),
                       c("Pvalue", "R2Adj.select")]
  r2 <- r2[!is.na(r2$R2Adj.select), ]
  write.csv(r2, paste0(output, "_selected-mems.csv"))
  # if best is mstree and second best is as good, take second best
  isBestMstree <- grepl("MST", rownames(r2)[1])
  isSecondBestAsGood <- (round(r2$R2Adj.select[1], 2) - round(r2$R2Adj.select[2], 2)) < 1e-2
  isSecondBestMstree <- grepl("MST", rownames(r2)[2])
  if (isBestMstree & isSecondBestAsGood & !isSecondBestMstree) {
    sel$best.id <- which(names(cands) == rownames(r2)[2])
    names(sel$best.id) <- rownames(r2)[2]
    sel$best <- mem.select(y, cands[[rownames(r2)[2]]])
    message("Took second best instead: ", names(sel$best.id))
  }
  # plot spatial graph and weighting function
  mode <- c("Delaunay", "Gabriel", "Relative", "MST", "Dnear0.3")
  mode <- switch (which(mode == strsplit(names(sel$best.id), "_")[[1]][1]),
                  "delauney", "gabriel", "relative", "mstree", "dnear")
  if (!is.na(mode)) {
    form <- c("Linear", "Up", "Down")
    par <- strsplit(names(sel$best.id), "_")[[1]][3]
    dmax <- round(max(unlist(nbdists(cands[[sel$best.id]]$neighbours, landscape))), 2)
    form <- switch (which(form == strsplit(names(sel$best.id), "_")[[1]][2]),
                    paste0("1 - x / ", dmax),
                    paste0("1 / (x ^ ", par, ")"),
                    paste0("1 - (x / ", dmax, ") ^ ", par))
    pdf(paste0(output, "_best-swm.pdf"), width = 8, height = 4)
    lw <- weight_ls(
      landscape,
      mode = mode,
      dist_formula = form,
      style = style,
      d2 = d2,
      plot = TRUE
    )
    dev.off()
  }
  # best MEMs summary
  write.csv(sel$best$summary, paste0(output, "_best-swm-summary.csv"))
  # best MEMs table
  write.csv(sel$best$MEM.select, paste0(output, "_best-swm-mems.csv"))
  # plot MEMs spatially
  pdf(paste0(output, "_plot-mems.pdf"), width = 4, height = 4)
  par(mfrow = c(2, 2))
  for (i in intersect(1:4, seq_len(ncol(sel$best$MEM.select)))) {
    s.value(landscape,
            sel$best$MEM.select[, i],
            method = "squaresize",
            pch = 20)
  }
  dev.off()
  # table of MEMs and pagerank centrality
  g <- graph_from_adjacency_matrix(
    nb2mat(cands[[names(sel$best.id)]]$neigh),
    weighted = TRUE
  )
  g <- as.undirected(g)
  un <- graph_from_adjacency_matrix(
    nb2mat(cands[[names(sel$best.id)]]$neigh) > 0
  )
  un <- as.undirected(un)
  centr <- as.data.frame(sel$best$MEM.select)
  # unweighted centralities
  centr$eigen <- eigen_centrality(un)$vector
  centr$degree <- centr_degree(un)$res
  centr$closeness <- centr_clo(un)$res
  centr$betweeness <- centr_betw(un)$res
  centr$eccentricity <- eccentricity(un)
  centr$pagerank <- page_rank(un)$vector
  centr$subgraph <- subgraph_centrality(un)
  # weighted centralities
  centr$w.eigen <- eigen_centrality(g)$vector
  centr$w.degree <- centr_degree(g)$res
  centr$w.closeness <- centr_clo(g)$res
  centr$w.betweeness <- centr_betw(g)$res
  centr$w.eccentricity <- eccentricity(g)
  centr$w.pagerank <- page_rank(g)$vector
  centr$w.subgraph <- subgraph_centrality(g)
  # save to file
  write.csv(centr, paste0(output, "_mems-centralities.csv"))
  return(list(
    mems = sel,
    centrality = centr,
    style = style
  ))
}
