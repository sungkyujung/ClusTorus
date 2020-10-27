#' K-Means Clustering to K-Spheres Clustering on Torus
#'
#' \code{kmeans.kspheres} prepares the parameters for conformity scores
#'   which are derived by k-means clustering on torus.
#'
#' @param data data n x d matrix of toroidal data on \eqn{[0, 2\pi)^2}
#' @param centers either the number of clusters or a set of initial
#'   cluster centers. If a number, a random set of row in x is
#'   chosen as the initial centers.
#' @param type character which must be "identical" or "various".
#'   If "identical", the radii of k-spheres are identical.
#'   If "various", the radii of k-spheres may be different.
#'
#' @return returns a \code{sphere.param} object,
#'   containing all values which determines the shape and
#'   location of spheres.
#'
#' @export
#' @seealso \code{\link[ClusTorus]{ehat.eval}},
#'   \code{\link[ClusTorus]{kmeans.torus}}
#' @references S. Jung, K. Park, and B. Kim (2020),
#'   "Clustering on the torus by conformal prediction", and
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
#' kmeans.kspheres(data, centers = 3, type = "various")
#' }
kmeans.kspheres <- function(data, centers = 10,
                            type = c("identical", "various")){

  # Returns a sphere.param object, containing all values which determines
  # the shape and location of spheres

  # type determines kmeans-fitting method. If "identical", the radii of
  # shperes are the same, and if not, the radii may be different.
  type <- match.arg(type)
  d <- ncol(data)
  n <- nrow(data)

  sphere.param <- list(mu1 = NULL, mu2 = NULL, Sigmainv = NULL, c = NULL)

  # Use extrinsic kmeans clustering for initial center points.
  # centers is given as a number, in default, but it may also be given
  # as a matrix which indicates the toroidal points.
  kmeans.out <- kmeans.torus(data, centers = centers,
                             method = "extrinsic")

  centers <- kmeans.out$centers
  J <- nrow(centers)
  # -------------- initializing ----------------

  # 1. identical spheres
  # In fact, the initialized parameters are for the identical case.
  sphere.param$mu1 <- centers[, 1]
  sphere.param$mu2 <- centers[, 2]
  sphere.param$c <- rep(0, J)

  for(j in 1:J){
    sphere.param$Sigmainv[[j]] <- diag(d)
  }


  # 2. various spheres
  if (type == "various"){

    for(j in 1:J){

      # if the size of cluster is 1, the cluster contains only one point.
      if (kmeans.out$size[j] == 0) { next }

      nj <- kmeans.out$size[j]
      pi_j <- nj / n
      sigma_j <- kmeans.out$withinss[j] / nj

      sphere.param$c[j] <- 2 * log(pi_j) - d * log(sigma_j)
      sphere.param$Sigmainv[[j]] <- diag(d) / sigma_j
    }
  }

  return(sphere.param)
}
