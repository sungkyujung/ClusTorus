#' Centers for explicit kmeans clustering on torus
#'
#' \code{centers.torus} returns a matrix, containing toroidal angles
#'   which are given by explicit kmeans algorithm on torus.
#'
#' @param data n x 2 matrix of toroidal data on \eqn{[0, 2\pi)^2}
#' @param k a scalar which determines the number of centers
#' @return a matrix containing toroidal angles
#' @export
#' @references S. Jung, K. Park, and B. Kim (2020),
#'   "Clustering on the torus by conformal prediction"
#' @examples
#' \dontrun{
#' ## mean vectors
#'
#' Mu1 <- c(3, 0)
#' Mu2 <- c(2, 2)
#' Mu3 <- c(1, 4)
#'
#' ## covariance matrices
#'
#' Sigma1 <- matrix(c(0.1, 0.05, 0.05, 0.2), 2, 2)
#' Sigma2 <- matrix(c(0.1, 0, 0, 0.01), 2, 2)
#' Sigma3 <- matrix(c(0.01, 0, 0, 0.1), 2, 2)
#'
#' ## 2-dimensional multivariate normal data wrapped with toroidal space
#' require(MASS)
#' data <- rbind(mvrnorm(n=70, Mu1, Sigma1),
#'               mvrnorm(n=50, Mu2, Sigma2),
#'               mvrnorm(n=50, Mu3, Sigma3))
#' data <- on.torus(data)
#'
#' centers.torus(data, k = 10)
#' }

centers.torus <- function(data, k = 10){
  kmeans.out <- ClusterR::KMeans_rcpp(cbind(cos(data),sin(data)), clusters = k)

  centroids <- kmeans.out$centroids
  centers <-cbind(atan2(centroids[, 3], centroids[, 1]), atan2(centroids[, 4], centroids[, 2]))

  centers <- on.torus(centers)
  return(centers)
}
