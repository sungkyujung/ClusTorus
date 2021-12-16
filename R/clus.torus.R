#' Clustering on the torus by conformal prediction
#'
#' \code{clus.torus} returns clustering results of data on the torus based on
#'   inductive conformal prediction set
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#'   or \eqn{[-\pi, \pi)^d}
#' @param split.id a n-dimensional vector consisting of values 1 (estimation)
#'   and 2(evaluation)
#' @param model A string. One of "mixture" and "kmeans" which
#'   determines the model or estimation methods. If "mixture", the model is based
#'   on the von Mises mixture, fitted
#'   with an EM algorithm. It supports the von Mises mixture and its variants
#'   based conformity scores. If "kmeans", the model is also based on the von
#'   Mises mixture, but the parameter estimation is implemented with the
#'   elliptical k-means algorithm illustrated in Appendix. It supports the
#'   log-max-mixture based conformity score only. If the
#'   dimension of data space is greater than 2, only "kmeans" is supported.
#'   Default is \code{method = "kmeans"}.
#' @param mixturefitmethod A string. One of "circular", "axis-aligned", and
#'   "general" which determines the constraint of the EM fitting. Default is
#'   "axis-aligned". This argument only works for \code{method = "mixture"}.
#' @param kmeansfitmethod A string. One of "general", ellipsoids",
#'   "heterogeneous-circular" or "homogeneous-circular". If "general", the
#'   elliptical k-means algorithm with no constraint is used. If "ellipsoids",
#'   only the one iteration of the algorithm is used. If"heterogeneous-circular",
#'   the same as above, but with the constraint that ellipsoids must be spheres.
#'   If "homogeneous-circular", the same as above but the radii of the spheres are
#'   identical. Default is "general". This argument only works for method = "kmeans".
#' @param J the number of components for mixture model fitting. If \code{J} is a vector,
#'   then \code{hyperparam.torus} is used to choose optimal \code{J}. If
#'  \code{J == NULL}, then \code{J = 4:30} is used.
#' @param level a scalar in \eqn{[0,1]}. The level of the conformal prediction set
#'  used for clustering. If \code{level == NULL}, then \code{hyperparam.alpha} is
#'  used to choose optimal \code{level}
#' @param option A string. One of "elbow", "risk", "AIC", or "BIC", which determines the
#'  criterion for the model selection. "risk" is based on the negative log-likelihood, "AIC"
#'  for the Akaike Information Criterion, and "BIC" for the Bayesian Information Criterion.
#'  "elbow" is based on minimizing the criterion used in Jung et. al.(2021).
#'  This argument is only used if \code{J} is a vector or \code{NULL}.
#' @param ... Further arguments that will be passed to \code{icp.torus} and
#'   \code{hyperparam.torus}
#' @return clustering assignment for data
#' \describe{
#'   \item{\code{cluster.id.by.log.density}}{cluster assignment result based on approximate log-density.}
#'   \item{\code{cluster.id.by.posterior}}{cluster assignment result based on the posterior probability.}
#'   \item{\code{cluster.id.outlier}}{cluster assignment result which regards data not included in conformal prediction set
#'   as outliers.}
#'   \item{\code{cluster.id.by.Mah.dist}}{cluster assignment result based on Mahalanobis distance.}
#'   \item{\code{level}}{used level which determines the size of clusters(conformal prediction set).}
#'   \item{\code{data}}{input data which are assigned to each cluster.}
#'   \item{\code{icp.torus}}{\code{icp.torus} object which is used for cluster assignment.}
#' }
#' @export
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction",
#' @examples
#' \donttest{
#' data <- toydata2[, 1:2]
#' n <- nrow(data)
#' clus.torus(data = data, model = "kmeans", kmeansfitmethod = "general", J = 5:30, option = "risk")
#' }

clus.torus <- function(data, split.id = NULL,
                       model = c("kmeans", "mixture"),
                       mixturefitmethod = c("axis-aligned","circular","general"),
                       kmeansfitmethod = c("general", "homogeneous-circular",
                                           "heterogeneous-circular",
                                           "ellipsoids"),
                       J = NULL,
                       level = NULL,
                       option = NULL,...){

  # this line prevents using "kde" for clustering.
  model <- match.arg(model)

  if(is.null(J)){J = 4:30}
  ind.J.select <- length(J) > 1

  if(is.null(level)){
    ind.alpha.select <- TRUE
  }else{
    ind.alpha.select <- FALSE
    if (length(level)>1){
      level = level[1]
      warning("Level must be a single number. Use the first element of argument level.")
    }
    if (level < 0 || level > 1) {
      level <- NULL
      ind.alpha.select <- TRUE
      warning("Level must be numeric and between 0 and 1. Reset as level = NULL (default)")
    }
  }

  model <- match.arg(model)
  mixturefitmethod <- match.arg(mixturefitmethod)
  kmeansfitmethod <- match.arg(kmeansfitmethod)


  if(ind.J.select){
    # J is a vector
    icp.torus.objects <- icp.torus(data, split.id = split.id,
                              model = model,
                              mixturefitmethod = mixturefitmethod,
                              kmeansfitmethod = kmeansfitmethod,
                              J = J,...)
    if (ind.alpha.select){
      # choose both J and alpha
      hyperparam.out <- hyperparam.torus(icp.torus.objects, option = option, ...)
      clust.out <- cluster.assign.torus(hyperparam.out)
      return(structure(list(cluster.obj = clust.out,
                  icp.torus = hyperparam.out$icp.torus,
                  hyperparam.select = hyperparam.out), class = "clus.torus"))
    }else{
      hyperparam.J.out <- hyperparam.J(icp.torus.objects, option = option)
      clust.out <- cluster.assign.torus(hyperparam.J.out$icp.torus, level = level)

      return(structure(list(cluster.obj = clust.out,
                  icp.torus = hyperparam.J.out$icp.torus,
                  hyperparam.select = hyperparam.J.out), class = "clus.torus"))
    }



  }else{

    # J is a scalar
    icp.torus <- icp.torus(data, split.id = split.id,
                                 model = model,
                                 mixturefitmethod = mixturefitmethod,
                                 kmeansfitmethod = kmeansfitmethod,
                                 J = J,...)


    if(ind.alpha.select){
      alpha.out <- hyperparam.alpha(icp.torus, ...)
      level = alpha.out$alphahat
    }else{
      alpha.out <- NULL
    }

    clust.out <- cluster.assign.torus(icp.torus, level = level)
    return(structure(list(cluster.obj = clust.out, icp.torus = icp.torus, hyperparam.select = alpha.out), class = "clus.torus"))
  }






}
