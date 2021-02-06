#' Transform the angular data to be on principal interval
#'
#' \code{on.torus} transforms d-dimensional angular data
#'   to be on \eqn{[0, 2\pi)^d}.
#'
#' @param x d-dimensional angular data(vector or matrix) whose unit is the radian.
#' @return d-dimensional radian-unit angular data on \eqn{[0, 2\pi)^d}.
#' @export
#' @examples
#' data <- SARS_CoV_2$tbl[1:200, 1:2]
#' data <- data * pi / 180
#'
#' on.torus(data)
on.torus <- function(x){
#  y <- x - 2 * pi * floor(x/2/pi + 1/2) # for [-pi, pi)^2
  y <- x - 2 * pi * floor(x/2/pi)
#  default version
#  y <- x
#  y[,1] <- x[,1] - 2*pi*floor(x[,1]/2/pi)
#  y[,2] <- x[,2] - 2*pi*floor(x[,2]/2/pi)
  return(y)
}
