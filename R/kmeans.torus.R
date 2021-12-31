#' K-Means Clustering on Torus
#'
#' \code{kmeans.torus} implements extrinsic k-means clustering
#'  on toroidal space.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#' @param centers either the number of clusters or a set of initial
#'   cluster centers. If a number, a random set of row in x is
#'   chosen as the initial centers.
#' @param ... Further arguments for \code{\link[stats]{kmeans}}.
#'
#' @details In Euclidean space, we know that the total sum of squares
#'   is equal to the summation of the within cluster sum of squares and
#'   the between cluster centers sum of squares. However, toroidal space
#'   does not satisfy the property; the equality does not hold. Thus,
#'   you need to be careful to use the sum of squares.
#'
#'   Extrinsic k-means algorithm uses the ambient space for \eqn{[0, 2\pi)^d}.
#'   Each datum is transformed to a vector in 2d-dimensional
#'   Euclidean space, whose elements are sine and cosine values of the datum,
#'   then a usual k-means algorithm is applied to transformed data.
#'
#' @return returns a \code{kmeans} object, which contains
#' \describe{
#'   \item{\code{extrinsic.results}}{extrinsic k-means clustering results using ordinary kmeans algorithm.}
#'   \item{\code{centers}}{A matrix of cluster centers.}
#'   \item{\code{membership}}{A vector of integers indicating the cluster to which each point is allocated.}
#'   \item{\code{size}}{The number of points in each cluster.}
#'   \item{\code{withinss}}{Vector of within-cluster sum of squares, one component per cluster.}
#'   \item{\code{totss}}{The total sum of squares, based on angular distance.}
#'   \item{\code{betweenss}}{The between-cluster sum of squares.}
#' }
#'
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
#'
#'   Gao, Y., Wang, S., Deng, M., & Xu, J. (2018). RaptorX-Angle: real-value prediction of protein backbone dihedral angles through a hybrid method of clustering and deep learning. \emph{BMC bioinformatics}, 19(4), 73-84.
#' @seealso \code{\link[stats]{kmeans}}, \code{\link[ClusTorus]{ang.minus}}, \code{\link[ClusTorus]{ang.dist}}
#' @export
#' @examples
#' data <- ILE[1:200, 1:2]
#'
#' kmeans.torus(data, centers = 2,
#'              iter.max = 100, nstart = 1)


kmeans.torus <- function(data, centers = 10, ...){

  # prepare for d - dimensional expansion
  n <- nrow(data)
  d <- ncol(data)

  kmeans <- list(data = data, centers = NULL,
                 membership = NULL, totss = NULL, withinss = NULL,
                 betweenss = NULL, size = NULL, extrinsic.results = NULL)

  # case for centers given as points on toroidal space
  if(length(centers) > 1){

    centroids <- cbind(cos(centers), sin(centers))
    kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = centroids, ...)

  }

  # case for centers given as a number or not given
  else if (length(centers) == 1){

    kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = centers, ...)
  }

  kmeans$extrinsic.results <- kmeans.out

  # calculate kmeans centers on torus
  centroids <- kmeans.out$centers
  centers <- atan2(centroids[, (d + 1):(2 * d)], centroids[, 1:d])
  centers <- on.torus(centers)

  kmeans$centers <- centers
  kmeans$membership <- kmeans.out$cluster
  kmeans$size <- kmeans.out$size

  J <- max(kmeans.out$cluster)

  # since the data are on T^d, totss is not the summation of
  # withinss and betweenss

  # -------- calculate withinss -----------

  # initialize with 0 vector
  kmeans$withinss <- rep(0, J)

  for(j in 1:J){

    # if the size of cluster is 1, withinss is 0
    if (kmeans$size[j] == 1) { next }

    # the case for the cluster size larger than 1
    nj <- length(data[kmeans$membership == j, ])
    j.mean <- wtd.stat.ang(data[kmeans$membership == j, ], rep(1, nj) / nj)$Mean


    kmeans$withinss[j] <- sum(tor.minus(data[kmeans$membership == j, ], j.mean)^2)
  }

  # --------- calculate totss ------------

  tot.mean <- wtd.stat.ang(data, rep(1, n) / n)$Mean
  kmeans$totss <- sum(tor.minus(data, tot.mean)^2)

  # --------- calculate betweenss ----------

  # if there is only one cluster, then center.distmat must be 0
  center.distmat <- ifelse(is.vector(centers), 0, ang.pdist(kmeans$centers))
  kmeans$betweenss <- sum(center.distmat^2)

  return(structure(kmeans, class = "kmeans.torus"))
}
