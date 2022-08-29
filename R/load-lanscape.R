#' @title Load Landscape
#'
#' @param ls_file file path of the food web adjacency table.
#'
#' @return a matrix with xy coordinates.
#'
#' @examples
#' \dontrun{
#' list.files("data/landscapes", full.names = TRUE)
#' i <- 1
#' path <- file.path("data", "landscapes")
#' ls_file <- file.path(path, paste0("RGG_", i, ".out"))
#' ls <- load_landscape(ls_file)
#' }
load_landscape <- function(ls_file) {
  # landscape type
  if (grepl("MI", ls_file)) {
    message("Mainland-island landscape")
  } else if (grepl("RGG", ls_file)) {
    message("Regular random landscape")
  } else {
    message("Small world landscape")
  }
  # load table
  patches <- read.csv(ls_file)
  # check area size is constant
  with_reporter(
    FailReporter,
    expect_equal(unique(patches$Area), 1e7)
  )
  xy <- patches[, c("X.Koord", "Y.Koord")]
  xy <- as.matrix(xy)
  colnames(xy) <- c("x", "y")
  return (xy)
}
