#' @title Default Color Palette
#'
#' @param n length of palette.
#'
#' @return String vector with colors.
webpal <- function(n) {
  cols <- c("#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF",
            "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850")
  ans <- colorRampPalette(cols)(n)
  return (ans)
}

#' @title Weighted Landscape
#'
#' @param xy landscape matrix.
#' @param mode method to construct the spatial network: "dnear", "knear",
#'   "delauney", "gabriel", "relative", "mstree".
#' @param d1 minimum distance for dnear.
#' @param d2 maximum distance for dnear.
#' @param k number of neighbors for knear.
#' @param dist_formula string of the function used to define the spatial
#'   weighting matrix (SWM).
#' @param style normalization style for SWM.
#' @param plot TRUE/FALSE.
#'
#' @return a list with
#'
#' @details Suggested dist_formulae:
#'   - "1 - (x / max(x))"
#'   - "1 - (x / max(x))^alpha", alpha in (2, 10)
#'   - "1 / x^beta", beta in (1, 10)
#'
#' @examples
#' \dontrun{
#' fw <- load_foodweb("data/foodwebs/web_1.out", "data/foodwebs/BodyMass_1.out")
#' ls <- load_landscape("data/landscapes/MI_1.out")
#' lsw <- weight_ls(ls, "dnear", d2 = max(V(fw)$dispersal), dist_formula = "1/x", plot = TRUE)
#' lsw <- weight_ls(ls, "delauney", plot = TRUE)
#' lsw <- weight_ls(ls, "knear", k = 1, plot = TRUE)
#' }
weight_ls <- function(
    xy,
    mode,
    d1 = 0,
    d2 = NULL,
    k = 0,
    dist_formula = "1 / x",
    style = "S",
    plot = FALSE
) {
  # some input checks
  if (!mode %in% c("dnear", "knear", "delauney", "gabriel", "relative", "mstree")) stop("mode not found")
  if (mode == "dnear" & (is.null(d1) | is.null(d2))) stop("dnear: specify d1 and d2")
  if (mode == "knear" & k == 0) stop("knear: specify k")
  nb <- switch(
    which(c("dnear", "knear", "delauney", "gabriel", "relative", "mstree") == mode),
    dnearneigh(xy, d1, d2),
    knn2nb(knearneigh(xy, k), sym = TRUE),
    tri2nb(xy),
    graph2nb(gabrielneigh(xy), sym = TRUE),
    graph2nb(relativeneigh(xy), sym = TRUE),
    mst.nb(dist(xy))
  )
  # spatial weighting matrix
  dlist <- nbdists(nb, xy)
  # glist <- lapply(dlist, function(x) 1 / x)
  glist <- lapply(dlist, function(x) eval(parse(text = dist_formula)))
  lw <- nb2listw(
    nb,
    glist,
    style,
    zero.policy = TRUE
  )
  swm <- listw2mat(lw)
  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    par(mar = c(4, 4, 0, 0))
    layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE),
           heights = c(.3, 1))
    plot.new()
    text(.5, .5, paste0("Neighbors: ", mode, "     ",
                        "Style: ", style, "     ",
                        "Distance formula :", dist_formula),
         cex = 1.5)
    plot(nb, xy,
         cex = 2.5,
         col = "steelblue", lw = 1.5)
    points(xy, pch = 20, cex = 3, col = "grey20")
    text(xy[, 1], xy[, 2], seq_len(nrow(xy)), col = "gold", cex = 1)
    x <- unlist(dlist)
    y <- unlist(lw$weights)
    fit <- loess(y ~ x)
    plot(x, y,
         pch = 20, cex = 1.5,
         col = adjustcolor("grey20", alpha.f = .5),
         frame = FALSE,
         xlab = "Distance", ylab = "Weight")
    lines(sort(x), fit$fitted[order(x)], lw = 2, col = "darkred")
    par(oldpar)
  }
  # return
  ans <- list(nb = nb, lw = lw, swm = swm)
  return (ans)
}
