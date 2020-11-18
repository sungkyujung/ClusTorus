#' Toroidal subtraction
#'
#' \code{tor.minus} computes angular subtraction bewtween n x d toroidal data and
#'   a d dimensional vector.
#'
#' @param data n x d matrix of toroidal data
#' @param mu a d-dimensinal vector
#' @return angular subtraction bewtween n x d toroidal data and
#'   a d dimensional vector.
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
  # data is n x d matrix of toroidal values, mu is a d-vector.
  # returns toroidal subtraction
  # cbind(
  #   ang.minus(data[,1], mu[1]),
  #   ang.minus(data[,2], mu[2]))
  
  # for higher dimension -----------
  t(apply(data, 1, ang.minus, Mu))
  
}
