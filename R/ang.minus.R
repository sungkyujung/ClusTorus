#' Angular subtraction
#'
#' \code{ang.minus} computes element-wise angular subtraction defined as
#' \deqn{x -  y := Arg(exp(i(x-y)))}
#'
#' @param x,y angular data(scalar or vector) whose elements are in \eqn{[0, 2\pi)}
#' @return returns a scalar or a vector whose elements are in
#'   \eqn{[-\pi, \pi)}.
#' @export
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. The Annals of Applied Statistics, 15(4), 1583-1603.
#' @examples
#' x <- c(pi/2, 0)
#' y <- c(pi, pi/3)
#'
#' ang.minus(x, y)
ang.minus <- function(x,y){
  # x "minus" y if x and y are [0, 2*pi)
  x <- on.torus(x)
  y <- on.torus(y)
  t <- rbind(x - y, x - y + 2 * pi, x - y - 2 * pi)

  tind <- apply(abs(t), 2, which.min)

  tt <- t[1, ]
  tt[tind == 2] <- t[2, tind == 2]
  tt[tind == 3] <- t[3, tind == 3]

  tt
}
