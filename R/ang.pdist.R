#' Pairwise L2 angular distance
#'
#' \code{ang.pdist} computes pairwise angular distances matrix.
#'
#' @param data n x d angular data on \eqn{[0, 2\pi)^d}
#' @return \code{ang.pdist} returns pairwise angular distances matrix
#'   with the class \code{dist}
#' @seealso \code{\link{ang.dist}}
#' @export
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#' @examples
#' data <- matrix(c(pi/3, pi/3, pi/2,
#'                  pi, pi/4, pi/2,
#'                  0, pi/3, pi/6),
#'                ncol = 3, byrow = TRUE)
#'
#' ang.pdist(data)
ang.pdist <- function(data){
  # assuming that data are n x d angular data on [0, 2pi)^d
  # computes L2 angular distance
  data <- on.torus(data)

  n <- nrow(data)
  d <- ncol(data)
  pdistmat <- matrix(0, ncol = n, nrow = n)

  for (i in 1:n-1){
    ad <- 0
    for (dd in 1:d){
      ad <- ad + (ang.dist(data[i, dd], data[(i + 1):n, dd]))^2
    }
    ad <- sqrt(ad)
    pdistmat[i, (i + 1):n] <- ad
    pdistmat[(i + 1):n, i] <- ad
  }
  stats::as.dist(pdistmat)
}
