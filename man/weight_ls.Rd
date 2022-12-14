% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weight-landscape.R
\name{weight_ls}
\alias{weight_ls}
\title{Weighted Landscape}
\usage{
weight_ls(
  xy,
  mode,
  d1 = 0,
  d2 = NULL,
  k = 0,
  dist_formula = "1 / x",
  style = "S",
  plot = FALSE
)
}
\arguments{
\item{xy}{landscape matrix.}

\item{mode}{method to construct the spatial network: "dnear", "knear",
"delauney", "gabriel", "relative", "mstree".}

\item{d1}{minimum distance for dnear.}

\item{d2}{maximum distance for dnear.}

\item{k}{number of neighbors for knear.}

\item{dist_formula}{string of the function used to define the spatial
weighting matrix (SWM).}

\item{style}{normalization style for SWM.}

\item{plot}{TRUE/FALSE.}
}
\value{
a list with
}
\description{
Weighted Landscape
}
\details{
Suggested dist_formulae:
  - "1 - (x / max(x))"
  - "1 - (x / max(x))^alpha", alpha in (2, 10)
  - "1 / x^beta", beta in (1, 10)
}
\examples{
\dontrun{
fw <- load_foodweb("data/foodwebs/web_1.out", "data/foodwebs/BodyMass_1.out")
ls <- load_landscape("data/landscapes/MI_1.out")
lsw <- weight_ls(ls, "dnear", d2 = max(V(fw)$dispersal), dist_formula = "1/x", plot = TRUE)
lsw <- weight_ls(ls, "delauney", plot = TRUE)
lsw <- weight_ls(ls, "knear", k = 1, plot = TRUE)
}
}
