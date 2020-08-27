#' Transform the angular values to be on principal interval
#'
#' \code{on.torus()} transforms data to be on \eqn{[0, 2\pi)^2}.
#'
#' @param x Matrix-formed numeric data which has 2 columns.
#' @export
#' @examples
#' \dontrun{
#' x <- matrix(c(10/3 * pi, 5/4 * pi), ncol = 2, byrow = TRUE)
#'
#' on.torus(x)
#' }
on.torus <- function(x){
#  y <- x - 2 * pi * floor(x/2/pi + 1/2) # for [-pi, pi)^2
  y <- x - 2 * pi * floor(x/2/pi)
#  default version
#  y <- x
#  y[,1] <- x[,1] - 2*pi*floor(x[,1]/2/pi)
#  y[,2] <- x[,2] - 2*pi*floor(x[,2]/2/pi)
  return(y)
}
