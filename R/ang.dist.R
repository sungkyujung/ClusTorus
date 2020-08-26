#' Angular distances
#'
#' \code{ang.dist()} computes angular distances.
#'
#' @param x,y angular data on \eqn{[-\pi, \pi)}
#' @export
#' @examples
#' \dontrun{
#' x <- pi/3
#' y <- -pi/3
#'
#' ang.dist(x, y)
#'}
ang.dist <- function(x,y){
  # dist(x,y) if x and y are [0, 2*pi]
  apply((rbind( abs(x - y),
                x + 2*pi - y,
                y + 2*pi - x)),2,min)
}
