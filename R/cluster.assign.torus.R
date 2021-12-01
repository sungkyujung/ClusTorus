#' Clustering by connected components of ellipsoids
#'
#' \code{cluster.assign.torus} returns clustering assignment for data
#'   given \code{icp.torus} objects, which can be constructed with
#'   \code{icp.torus.score}.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}.
#' @param icp.torus an object containing all values to compute the conformity
#'   score, which will be constructed with \code{icp.torus.score}.
#' @param level a scalar in \eqn{[0,1]}. Default value
#'   is 0.1.
#' @return clustering assignment for data, given \code{icp.torus} objects
#' @export
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction",
#'
#'   I. Gilitschenski and U. D. Hanebeck,
#'   "A robust computational test for overlap of
#'   two arbitrary-dimensional ellipsoids in fault-detection of Kalman filters"
#' @seealso \code{\link[ClusTorus]{icp.torus.score}}
#' @examples
#' data <- toydata1[, 1:2]
#' icp.torus <- icp.torus.score(data, method = "kmeans",
#'                              kmeansfitmethod = "general",
#'                              J = 4, concentration = 25)
#' level <- 0.1
#'
#' cluster.assign.torus(data, icp.torus, level)
cluster.assign.torus <- function(data, icp.torus, level = 0.1){
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
  data <- on.torus(data)

  n2 <- icp.torus$n2
  ialpha <- ifelse((n2 + 1) * level < 1, 1, floor((n2 + 1) * level))

  cluster.obj <- list()
  # For kmeans to kspheres --------------------------------------------------

  if(icp.torus$method == "kde"){
    stop("method kde does not support clustering.")
  } else {
    t <- icp.torus$score_ellipse[ialpha]
    cluster.obj <- conn.comp.ellipse(icp.torus$ellipsefit, t)
    K <- cluster.obj$ncluster
    ehatj <- ehat.eval(data, icp.torus$ellipsefit)

    ehatj[,cluster.obj$componentid == 0] <- -Inf

    # maxj.id <- apply(ehatj, 1, which.max)
    maxj.id <- max.col(ehatj, ties.method = "first")
    cluster.id1 <- cluster.obj$componentid[maxj.id]
    cluster.obj$cluster.id.by.ehat <- cluster.id1

    partsum <- matrix(0, nrow = nrow(data), ncol = K)
    for(k in 1:K){
      ifelse(sum( cluster.obj$componentid == k) > 1,
             partsum[,k] <- rowSums(ehatj[,cluster.obj$componentid == k]),
             partsum[,k] <- ehatj[,cluster.obj$componentid == k])
    }
    # cluster.obj$mixture$cluster.id.by.partialsum <- apply(partsum, 1, which.max)
    cluster.obj$cluster.id.by.partialsum <- max.col(partsum, ties.method = "first")

    # ehat <- apply(ehatj,1,max)
    ehat <- do.call(pmax, as.data.frame(ehatj))

    cluster.id1[!(ehat >= icp.torus$score_ellipse[ialpha])] <- K+1
    cluster.obj$cluster.id.outlier <- cluster.id1

    # Use mahalanobis distance

    mah <- mah.dist2.eval(data, icp.torus$ellipsefit)
    mah[,cluster.obj$componentid == 0] <- Inf

    # maxj.id <- apply(mah, 1, which.min)
    maxj.id <- max.col(-mah, ties.method = "first")
    cluster.id1 <- cluster.obj$componentid[maxj.id]
    cluster.obj$cluster.id.by.Mah.dist <- cluster.id1

  }
  cluster.obj$level <- level
  cluster.obj$data <- data
  cluster.obj$icp.torus <- icp.torus
  return(structure(cluster.obj, class = "cluster.obj"))
}
