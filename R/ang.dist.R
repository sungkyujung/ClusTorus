#' Angular distance
#'
#' \code{ang.dist} computes element-wise angular distance between
#'  two angular values in \eqn{[0,2\pi)}.
#'
#' @param x,y angular data(both scalar or vector) whose elements are in \eqn{[0, 2\pi)}
#' @return angular data (scalar or vector) whose elements are in \eqn{[0, 2\pi)}
#' @export
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
#' @examples
#' x <- c(pi/3, 0)
#' y <- c(pi/4, pi/2)
#'
#' ang.dist(x, y)
ang.dist <- function(x,y){
  # dist(x,y) if x and y are [0, 2*pi]
  apply((rbind( abs(x - y),
                x + 2*pi - y,
                y + 2*pi - x)),2,min)
}
