#' Kernel density estimation using circular von Mises distribution
#'
#' \code{kde.torus} returns a kde using independent multivariate von mises kernel.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#' @param eval.point N x N numeric matrix on \eqn{[0, 2\pi)^d}. Default input is
#'  \code{NULL}, which represents the fine grid points on \eqn{[0, 2\pi)^d}.
#' @param concentration positive number which has the role of \eqn{\kappa} of
#'   von Mises distribution. Default value is 25.
#' @return \code{kde.torus} returns N-dimensional vector of kdes evaluated at eval.point
#' @export
#' @seealso \code{\link{grid.torus}}
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#' @examples
#' data <- ILE[1:200, 1:2]
#'
#' kde.torus(data)

kde.torus <- function(data, eval.point = NULL,
                      concentration = 25){
  # Computes kde using independent bivariate von mises kernel.
  # evaluated at eval.point, a N x 2 matrix of points on torus
  # returns N-vector of kdes evaluated at eval.point
  # used in cp.torus.kde()
  if (!is.matrix(data)) {data <- as.matrix(data)}

  d <- ncol(data)
  n <- nrow(data)

  if (is.null(eval.point)){
    grid.size <- ifelse(d == 2, 100, floor(10^(6/d)))
    eval.point <- grid.torus(d = d, grid.size = grid.size)
  }

  N <- nrow(eval.point)

  summand <- rep(0, N)
  for (i in 1:n){
    summand <- summand + exp(concentration *
                               rowSums(cos(sweep(eval.point, 2, data[i, ], "-"))))
  }
  return(summand / n / (2 * pi * besselI(concentration, 0))^d)
}
