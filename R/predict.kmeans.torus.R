#' Prediction for Extrinsic Kmeans Clustering
#'
#' \code{pred.kmeans.torus} predicts the cluster for each
#'   data point.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}.
#' @param kmeans a \code{kmeans} object, which contains
#'   input data, cluster centers on torus, membership,
#'   total sum of squares, within cluster sum of squares,
#'   between cluster centers sum of squares, and the size of
#'   each cluster. See  \code{\link[ClusTorus]{kmeans.torus}}
#' @return a vector whose elements indicate the labels of predicted clusters.
#' @export
#' @references 'S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#' @seealso \code{\link[ClusTorus]{kmeans.torus}}
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
#' split.id <- sample(1:2, nrow(data), replace = TRUE)
#' data.train <- data[split.id == 1, ]
#' data.test <- data[split.id == 2, ]
#'
#' kmeans <- kmeans.torus(data.train, centers = 3,
#'              iter.max = 100, nstart = 1)
#'
#' pred.kmeans.torus(data.test, kmeans)
#' }
pred.kmeans.torus <- function(data, kmeans){
  data <- on.torus(data)

  extrinsic.results <- kmeans$extrinsic.results
  extrinsic.data <- cbind(cos(data), sin(data))

  pred.kmeans <- apply(extrinsic.data, 1, function(r)
    {which.min(colSums((t(extrinsic.results$centers) - r)^2))})

  return(pred.kmeans)
}
