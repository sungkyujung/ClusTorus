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
#' @references S. Jung, K. Park, and B. Kim (2020),
#'   "Clustering on the torus by conformal prediction"
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
  if(is.vector(data) || is.data.frame(data)){ data <- as.matrix(data) }
  if(dim(data)[2] == 1) {data <- t(data)}
  if(!is.vector(mu) && ncol(mu) == 1) {mu <- as.vector(mu)}
  if(!is.numeric(data)) {stop("Invaild data : data must be numeric.")}
  if(!is.numeric(mu)) {stop("Invaild data : mu must be numeric.")}
  # cbind(
  #  ang.minus(data[,1], mu[1]),
  #  ang.minus(data[,2], mu[2]))

  # for higher dimension -----------
  tor <- ang.minus(data[,1], mu[1])
  d <- length(mu)
  if (d > 1){
    for (i in 2:d){
      tor <- cbind(tor, ang.minus(data[,i], mu[i]))
    }
  }
  tor
}
