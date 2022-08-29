#' @title Load Time Series
#'
#' @param ts_path path of the time series
#'
#' @examples
#' \dontrun{
#' f <- "~/Proj/remilio/data/output/remilio/6/timeseries_MI_1.out"
#' }
load_timeseries <- function(ts_path) {
  nlines <- system(paste("wc -l", ts_path), intern = TRUE)
  nlines <- strsplit(nlines, " ")[[1]][1]
  nlines <- as.numeric(nlines)
  ts <- read.table(ts_path, skip = nlines - 1)
  time <- as.numeric(ts[1])
  biomass <- as.numeric(ts[-1])
  ts_names <- read.table(ts_path, nrows = 1)
  ts_names <- ts_names[-1]
  patch <- sapply(ts_names, function(x) strsplit(x, "P")[[1]][[2]])
  patch <- as.numeric(unlist(patch))
  sp <- sapply(ts_names, function(x) strsplit(x, "P")[[1]][[1]])
  sp <- as.numeric(gsub("B", "", unlist(sp)))
  # re-adjust for R indexing
  if (min(patch) == 0) patch <- patch + 1
  if (min(sp) == 0) sp <- sp + 1
  ans <- data.frame(
    patch = patch,
    species = sp,
    time = time,
    biomass = biomass
  )
  return (ans)
}

#' @title Local Community Composition
#'
#' @param timeseries data.frame created with load_ts().
#' @param foodweb igraph object created with load_foodweb().
#'
#' @return an igraph object.
local_webs <- function(timeseries, foodweb) {
  patch <- as.list(rep(NA, max(timeseries$patch)))
  for (i in seq_along(patch)) {
    #patch[[i]]
    biomass <- timeseries$biomass[timeseries$patch == i]
    sp <- which(biomass > 0)
    local_fw <- igraph::induced_subgraph(foodweb, sp)
    patch[[i]] <- local_fw
  }
  props <- lapply(patch, network)
  props <- do.call("rbind", props)
  props$patch <- seq_len(nrow(props))
  jaccard <- matrix(NA, nrow = length(patch), ncol = length(patch))
  for (i in seq_len(nrow(jaccard))) {
    for (j in seq_len(ncol(jaccard))) {
      species <- lapply(patch[c(i, j)], function(x) V(x))
      cup <- length(unique(unlist(species)))
      cap <- length(intersect(species[[1]], species[[2]]))
      jaccard[i, j] <- cap / cup
    }
  }
  ans <- list(
    webs = patch,
    jaccard = jaccard,
    properties = props
  )
  return (ans)
}
