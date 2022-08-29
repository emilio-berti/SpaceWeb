#' @title Weighted Landscape
#'
#' @param xy landscape matrix.
#' @param mode method to construct the spatial network: "dnear", "knear",
#'   "delauney", "gabriel", "relative".
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
  if (!mode %in% c("dnear", "knear", "delauney", "gabriel", "relative")) stop("mode not found")
  if (mode == "dnear" & (is.null(d1) | is.null(d2))) stop("dnear: specify d1 and d2")
  if (mode == "knear" & k == 0) stop("knear: specify k")
  nb <- switch(
    which(c("dnear", "knear", "delauney",
            "gabriel", "relative", "mstree") == mode),
    dnearneigh(xy, d1, d2),
    knn2nb(knearneigh(xy, k)),
    tri2nb(xy),
    graph2nb(gabrielneigh(xy)),
    graph2nb(relativeneigh(xy)),
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
         cex = 1,
         col = "steelblue", lw = 1.5)
    points(xy, pch = 20, col = "grey20")
    x <- unlist(dlist)
    y <- unlist(swm[swm > 0])
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
