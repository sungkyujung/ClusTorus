#' Clustering by connected components of ellipses
#'
#' \code{cluster.assign.torus} returns clustering assignment for data
#'   given \code{icp.torus} objects, which can be constructed with
#'   \code{icp.torus.score}.
#'
#' @param data n x 2 matrix of toroidal data on \eqn{[0, 2\pi)^2}.
#' @param icp.torus an object containing all values to compute the conformity
#'   score, which will be constructed with \code{icp.torus.score}.
#' @param level either a scalar or a vector, or even \code{NULL}. Default value
#'   is 0.1.
#' @return clustering assignment for data, given icp.torus objects
#' @export
#' @references 'S. Jung, K. Park, and B. Kim (2020),
#'   "Clustering on the torus by conformal prediction"
#' @seealso \code{\link[ClusTorus]{icp.torus.score}}
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
#' icp.torus <- icp.torus.score(data, method = "all",
#'                              mixturefitmethod = "general",
#'                              param = list(J = 4, concentration = 25))
#' level <- c(0.1, 0.08)
#'
#' cluster.assign.torus(data, icp.torus, level)
#' }
cluster.assign.torus <- function(data, icp.torus, level){
  # clustering by connected components of ellipses
  #
  # return clustering assignment for data, given icp.torus objects.
  # clustering assignment is given by
  # 1) max_j phatj (phatj = pi_j f_j(x)) for max-mixture (not implemented yet)
  # 2) max_k sum_{j in cluster k} phatj for max-mixture  (not implemented yet)
  # 3) max_j ehatj for ellipse-approx
  # 4) max_k sum_{j in cluster k} ehatj for ellipse-approx
  # two options
  # one: every point is assigned to clusters 1:K by either (1)-(4) above
  # two: outliers are assigned to cluster K+1.
  n2 <- icp.torus$n2
  ialpha <- floor( (n2 + 1) * level)

  # For max-ellipse ---------------------------------------------------------

  t <- icp.torus$mixture$score_ellipse[ialpha]
  cluster.obj <- conn.comp.ellipse(icp.torus$mixture$ellipsefit, t)
  K <- cluster.obj$ncluster
  ehatj <- ehat.eval(data, icp.torus$mixture$ellipsefit)

  ehatj[,cluster.obj$componentid == 0] <- -Inf

  maxj.id <- apply(ehatj, 1, which.max)
  cluster.id1 <- cluster.obj$componentid[maxj.id]
  cluster.obj$cluster.id.by.ehat <- cluster.id1

  partsum <- matrix(0, nrow = nrow(data), ncol = K)
  for(k in 1:K){
    ifelse(sum( cluster.obj$componentid == k) > 1,
           partsum[,k] <- rowSums(ehatj[,cluster.obj$componentid == k]),
           partsum[,k] <- ehatj[,cluster.obj$componentid == k])
  }
  cluster.obj$cluster.id.by.partialsum <- apply(partsum, 1, which.max)

  ehat <- apply(ehatj,1,max)

  cluster.id1[!(ehat >= icp.torus$mixture$score_ellipse[ialpha])] <- K+1
  cluster.obj$cluster.id.outlier <- cluster.id1

  # Use mahalanobis distance

  mah <- mah.dist2.eval(data, icp.torus$mixture$ellipsefit)
  mah[,cluster.obj$componentid == 0] <- Inf

  maxj.id <- apply(mah, 1, which.min)
  cluster.id1 <- cluster.obj$componentid[maxj.id]
  cluster.obj$cluster.id.by.Mah.dist <- cluster.id1


  return(cluster.obj)
}
