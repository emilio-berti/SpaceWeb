#' @title Load Food Web
#'
#' @param fw_file file path of the food web adjacency table.
#' @param sp_file file path of the species data.frame.
#'
#' @return an igraph object.
#'
#' @examples
#' \dontrun{
#' i <- 1
#' path <- file.path("data", "foodwebs")
#' sp_file <- file.path(path, paste0("BodyMass_", i, ".out"))
#' fw_file <- file.path(path, paste0("web_", i, ".out"))
#' fw <- load_foodweb(fw_file, sp_file)
#' }
load_foodweb <- function(fw_file, sp_file) {
  sp <- read.csv(sp_file)
  # craete graph
  adj <- read.table(fw_file)
  adj <- as.matrix(adj)
  adj <- t(adj)
  g <- graph_from_adjacency_matrix(adj)
  # species traits
  V(g)$mass <- sp$body.mass
  V(g)$dispersal <- sp$D_Max_i
  # hald saturation densities
  V(g)$N1 <- sp$Nutr.halfsat.dens1
  V(g)$N2 <- sp$Nurt.halfsat.dens2
  # role of species
  deg_out <- degree(g, mode = "out")
  deg_in <- degree(g, mode = "in")
  V(g)$role <- NA
  V(g)$role[which(deg_out > 0 & deg_in == 0)] <- "basal"
  V(g)$role[which(deg_out > 0 & deg_in > 0)] <- "consumer"
  V(g)$role[which(deg_out == 0 & deg_in > 0)] <- "top"
  # visualization pars
  V(g)$color <- NA
  V(g)$color[V(g)$role == "basal"] <- "chartreuse3"
  V(g)$color[V(g)$role == "consumer"] <- "goldenrod2"
  V(g)$color[V(g)$role == "top"] <- "tomato"
  #some checks
  test_basal <- identical(as.logical(sp$if.basal.spp), V(g)$role == "basal")
  with_reporter(
    FailReporter,
    expect_identical(as.logical(sp$if.basal.spp), V(g)$role == "basal")
  )
  with_reporter(
    FailReporter,
    expect_identical(sp$body.mass, V(g)$mass)
  )
  with_reporter(
    FailReporter,
    expect_identical(sp$D_Max_i, V(g)$dispersal)
  )
  return (g)
}
