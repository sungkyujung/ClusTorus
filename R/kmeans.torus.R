#' K-Means Clustering on Torus
#'
#' \code{kmeans.torus} implements extrinsic k-means clustering
#'  on toroidal space.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#' @param centers either the number of clusters or a set of initial
#'   cluster centers. If a number, a random set of row in x is
#'   chosen as the initial centers.
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
#' @references 'S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#'
#' @seealso \code{\link[stats]{kmeans}}, \code{\link[ClusTorus]{ang.minus}}
#' @export
#' @examples
#' data <- ILE[1:200, 1:2]
#'
#' kmeans.torus(data, centers = 2,
#'              iter.max = 100, nstart = 1)


kmeans.torus <- function(data, centers = 10,
                         iter.max = 100, nstart = 1){

  # prepare for d - dimensional expansion
  n <- nrow(data)
  d <- ncol(data)

  kmeans <- list(data = data, centers = NULL,
                 membership = NULL, totss = NULL, withinss = NULL,
                 betweenss = NULL, size = NULL, extrinsic.results = NULL)

  # case for centers given as points on toroidal space
  if(length(centers) > 1){

    centroids <- cbind(cos(centers), sin(centers))
    kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = centroids,
                         iter.max = iter.max, nstart = nstart)

  }

  # case for centers given as a number or not given
  else if (length(centers) == 1){

    kmeans.out <- kmeans(cbind(cos(data),sin(data)), centers = centers,
                         iter.max = iter.max, nstart = nstart)
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

  return(kmeans)
}
