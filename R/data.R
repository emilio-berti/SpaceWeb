#' Landscape
#'
#' xy coordinates for 60 patches for the Mainland-Island landscape
#'
#' @format A matrix with 60 rows and 2 columns:
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordiante}
#'   ...
#' }
"landscape"

#' Meta-foodweb
#'
#' A graph of the meta-foodweb
#'
#' @format An igraph object with 40 vertices and 218 edges.
"foodweb"

#' Time series
#'
#' One time series of species biomasses.
#'
#' @format A data frame with 2400 rows and 4 features:
#' \describe{
#'   \item{patch}{patch ID}
#'   \item{species}{species ID}
#'   \item{time}{integration time step}
#'   \item{biomass}{biomass of species}
#'   ...
#' }
"timeseries"
