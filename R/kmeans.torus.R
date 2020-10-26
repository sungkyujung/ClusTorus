#' K-Means Clustering on Torus
#'
#' \code{kmeans.torus} implements k-means clustering on toroidal
#'  space.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#' @param centers either the number of clusters or a set of initial
#'   cluster centers. If a number, a random set of row in x is
#'   chosen as the initial centers.
#' @param method a string one of "intrinsic" or "extrinsic".
#'   If "extrinsic", then extrinsic kmeans clustering method will be
#'   implemented. ("intrinsic" is not yet supported)
#' @param iter.max the maximum number of iterations
#' @param nstart if \code{centers} is a number, how many random sets
#'   should be chosen?
#'
#' @details In Euclidean space, we know that the total sum of squares
#'   is equal to the summation of the within cluster sum of squares and
#'   the between cluster centers sum of squares. However, toroidal space
#'   does not satisfy the property; the equality does not hold. Thus,
#'   you need to be careful to use the sum of squares.
#'
#' @return returns a \code{kmeans} object, which contains
#'   input data, cluster centers on torus, membership,
#'   total sum of squares, within cluster sum of squares,
#'   between cluster centers sum of squares, and the size of
#'   each cluster.
#'
#' @references 'S. Jung, K. Park, and B. Kim (2020),
#'   "Clustering on the torus by conformal prediction", and
#'   Jaehyeok Shin, Alessandro Rinaldo and Larry Wasserman (2019),
#'   "Predictive Clustering"
#'
#' @seealso \code{\link[stats]{kmeans}}, \code{\link[ClusTorus]{ang.minus}}
#' @export
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
#' kmeans.torus(data, centers = 3, method = "extrinsic",
#'              iter.max = 100, nstart = 1)
#' }


kmeans.torus <- function(data, centers = 10,
                         method = c("intrinsic", "extrinsic"),
                         iter.max = 100, nstart = 1){

  # prepare for d - dimensional expansion
  d <- ncol(data)
  method <- match.arg(method)

  kmeans <- list(data = data, centers = NULL,
                 membership = NULL, totss = NULL, withinss = NULL,
                 betweenss = NULL, size = NULL)

  # 1. Extrinsic kmeans
  if(method == "extrinsic"){

    # case for centers given as points on toroidal space
    if(length(centers) > 1){

      centers <- cbind(cos(centers), sin(centers))
      kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = centers,
                           iter.max = iter.max, nstart = nstart)

    }

    # case for centers given as a number or not given
    else if (length(centers) == 1){

      kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = centers,
                           iter.max = iter.max, nstart = nstart)
    }

    # calculate kmeans centers on torus
    centroids <- kmeans.out$centers
    centers <- atan2(centroids[, (d + 1):(2 * d)], centroids[, 1:d])
    centers <- on.torus(centers)

    kmeans$centers <- centers
    kmeans$membership <- kmeans.out$cluster
    kmeans$size <- kmeans.out$size

    k <- nrow(centers)

    # since the data are on T^d, totss is not the summation of
    # withinss and betweenss

    # calculate withinss
    for(j in 1:k){

      cluster.dev <- t(apply(data[kmeans$membership == j, ], 1, ang.minus,
                           y = colMeans(data[kmeans$membership == j, ])))

      kmeans$withinss[j] <- sum(cluster.dev^2)
    }

    # calculate totss
    tot.mean <- colMeans(data)
    tot.dev <- t(apply(data, 1, ang.minus, y = tot.mean))
    kmeans$totss <- sum(tot.dev^2)

    # calculate betweenss
    center.distmat <- ang.pdist(kmeans$centers)
    kmeans$betweenss <- sum(center.distmat^2)

    # 2. Intrinsic kmeans(not yet bulit)
  } else { stop("Intrinsic not yet implemented") }

  return(kmeans)
}