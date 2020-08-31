#' Toroidal subtraction
#'
#' \code{tor.minus} computes angular subtraction bewtween n x 2 toroidal data and
#'   a 2 dimensional vector.
#'
#' @param data n x 2 matrix of toroidal data
#' @param mu a 2-dimensinal vector
#' @return angular subtraction bewtween n x 2 toroidal data and
#'   a 2 dimensional vector.
#' @seealso \code{\link{ang.minus}}
#' @export
#' @examples
#' \dontrun{
#' data <- matrix(c(pi/3, pi/4, pi/2, pi/2),
#'                ncol = 2, byrow = TRUE)
#' mu <- c(pi/2, pi/2)
#'
#' tor.minus(data, mu)
#' }

tor.minus <- function(data, mu){
  # data is n x 2 matrix of toroidal values, mu is a 2-vector.
  # returns toroidal subtraction
  cbind(
    ang.minus(data[,1], mu[1]),
    ang.minus(data[,2], mu[2]))
}
