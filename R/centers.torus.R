#' Centers from extrinsic kmeans clustering on torus (depreciated)
#'
#' \code{centers.torus} returns a matrix, containing toroidal angles
#'   which are given by extrinsic kmeans algorithm on torus.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#' @param k a scalar which determines the number of centers
#' @return a matrix containing toroidal angles which indicates
#'   k centers.
#' @export
#' @references S. Jung, K. Park, and B. Kim (2020),
#'   "Clustering on the torus by conformal prediction" and
#'   Jaehyeok Shin, Alessandro Rinaldo and Larry Wasserman (2019),
#'   "Predictive Clustering"
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
  # kmeans.out <- ClusterR::KMeans_rcpp(cbind(cos(data),sin(data)), clusters = k)
  kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = k)
  centroids <- kmeans.out$centers
  centers <-atan2(centroids[, (d + 1):(2 * d)], centroids[, 1:d])

  centers <- on.torus(centers)
  return(centers)
}
