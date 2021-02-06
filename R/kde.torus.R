#' Kernel density estimation using circular von Mises distribution
#'
#' \code{kde.torus} returns a kde using independent bivariate von mises kernel.
#'
#' @param data n x 2 matrix of toroidal data on \eqn{[0, 2\pi)^2}
#' @param eval.point N x N numeric matrix on \eqn{[0, 2\pi)^2}. Default input is
#'  \code{grid.torus}.
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

kde.torus <- function(data, eval.point = grid.torus(),
                      concentration = 25){
  # Computes kde using independent bivariate von mises kernel.
  # evaluated at eval.point, a N x 2 matrix of points on torus
  # returns N-vector of kdes evaluated at eval.point
  # used in cp.torus.kde()
  if (!is.matrix(data)) {data <- as.matrix(data)}

  n <- nrow(data)
  N <- nrow(eval.point)

  summand <- rep(0, N)
  for (i in 1:n){
    summand <- summand + exp(concentration *
                               rowSums(cos(sweep(eval.point, 2, data[i, ], "-"))))
  }
  return(summand / n / (2 * pi * besselI(concentration, 0))^2)
}
