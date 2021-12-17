#' Selecting optimal level based on the runs of the number of clusters
#'
#' \code{hyperparam.alpha} evaluates the numbers of clusters for various
#'   levels, and select the optimal level based on the runs of the cluster numbers.
#'
#' @param icp.torus an object containing all values to compute the conformity
#'   score, which will be constructed with \code{icp.torus.score}.
#' @param alphavec either a scalar or a vector, or even \code{NULL} for the levels.
#'   Default value is \code{NULL}. If \code{NULL}, then \code{alphavec} is
#'   automatically generated as a sequence from 0 to \code{alpha.lim}.
#' @param alpha.lim a positive number lower than 1, which is the upper bound of
#'   Default is 0.15.
#' @return returns a \code{hyperparam.alpha} object which contains a \code{data.frame} for
#'   the numbers of clusters corresponding to the levels and the optimal level.
#' @export
#' @seealso \code{\link[ClusTorus]{hyperparam.J}}, \code{\link[ClusTorus]{hyperparam.torus}}
#'  \code{\link[ClusTorus]{icp.torus}}
#' @examples
#' \donttest{
#' data <- toydata2[, 1:2]
#' n <- nrow(data)
#' split.id <- rep(2, n)
#' split.id[sample(n, floor(n/2))] <- 1
#' icp.torus <- icp.torus(data, split.id = split.id, model = "kmeans",
#'                        kmeansfitmethod = "ge", init = "h",
#'                        J = 25, verbose = TRUE)
#' hyperparam.alpha(icp.torus)
#' }
hyperparam.alpha <- function(icp.torus, alphavec = NULL, alpha.lim = 0.15){
  if(is.null(icp.torus)) {stop("icp.torus object must be input.")}

  if(icp.torus$model == "mixture") {model <- "mixture"}
  else if(icp.torus$model == "kmeans") {model <- "kmeans"}
  else {stop("model kde is not supported.")}
  n2 <- icp.torus$n2

  if (is.null(alphavec) && alpha.lim > 1) {stop("alpha.lim must be less than 1.")}

  output <- list()
  out <- data.frame()

  if (is.null(alphavec)) {alphavec <- 1:floor(min(n2, 1000) * alpha.lim) / n2}

  n.alphavec <- length(alphavec)

  # 1. kmeans -----------------------------------------------------
  if (model == "kmeans"){

    for (ii in 1:n.alphavec){
      alpha <- alphavec[ii]
      ialpha <- ifelse((n2 + 1) * alpha < 1, 1, floor((n2 + 1) * alpha))
      t <- icp.torus$score_ellipse[ialpha]
      ncluster <- conn.comp.ellipse(icp.torus$ellipsefit, t)$ncluster

      out <- rbind(out, data.frame(alpha = alpha, ncluster = ncluster))
      if(ii%%10 == 0) cat(".")
    }
    cat("\n")

    nclusters.length <- rle(out$ncluster)$lengths
    length <- max(nclusters.length)
    length.index <- which.max(nclusters.length)
    length.sum <- ifelse(length.index == 1, 0, sum(nclusters.length[1:(length.index - 1)]))

    term.alpha <- out$alpha[(length.sum + 1):(length.sum + length)]
    alphahat <- stats::median(term.alpha)

    output$alpha.results <- out
    output$alphahat <- alphahat
  }

  # 2. mixture ----------------------------------------------------

  else if (model == "mixture"){

    for (ii in 1:n.alphavec){
      alpha <- alphavec[ii]
      ialpha <- ifelse((n2 + 1) * alpha < 1, 1, floor((n2 + 1) * alpha))
      t <- icp.torus$score_ellipse[ialpha]
      ncluster <- conn.comp.ellipse(icp.torus$ellipsefit, t)$ncluster

      out <- rbind(out, data.frame(alpha = alpha, ncluster = ncluster))
      if(ii%%10 == 0) cat(".")
    }
    cat("\n")

    nclusters.length <- rle(out$ncluster)$lengths
    length <- max(nclusters.length)
    length.index <- which.max(nclusters.length)
    length.sum <- ifelse(length.index == 1, 0, sum(nclusters.length[1:(length.index - 1)]))

    term.alpha <- out$alpha[(length.sum + 1):(length.sum + length)]
    alphahat <- stats::median(term.alpha)

    output$alpha.results <- out
    output$alphahat <- alphahat
  }
  return(structure(output, class = "hyperparam.alpha"))
}
