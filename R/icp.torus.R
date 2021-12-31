#' Conformity score for inductive prediction sets
#'
#' \code{icp.torus} prepares all values
#'   for computing the conformity score for specified methods.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#'   or \eqn{[-\pi, \pi)^d}
#' @param split.id a n-dimensional vector consisting of values 1 (estimation)
#'   and 2(evaluation)
#' @param model A string. One of "kde", "mixture", and "kmeans" which
#'   determines the model or estimation methods. If "kde", the model is based
#'   on the kernel density estimates. It supports the kde-based conformity score
#'   only. If "mixture", the model is based on the von Mises mixture, fitted
#'   with an EM algorithm. It supports the von Mises mixture and its variants
#'   based conformity scores. If "kmeans", the model is also based on the von
#'   Mises mixture, but the parameter estimation is implemented with the
#'   elliptical k-means algorithm illustrated in Appendix. It supports the
#'   log-max-mixture based conformity score only. If the
#'   dimension of data space is greater than 2, only "kmeans" is supported.
#'   Default is \code{model = "kmeans"}.
#' @param mixturefitmethod A string. One of "circular", "axis-aligned", and
#'   "general" which determines the constraint of the EM fitting. Default is
#'   "axis-aligned". This argument only works for \code{model = "mixture"}.
#' @param kmeansfitmethod A string. One of "general", ellipsoids",
#'   "heterogeneous-circular" or "homogeneous-circular". If "general", the
#'   elliptical k-means algorithm with no constraint is used. If "ellipsoids",
#'   only the one iteration of the algorithm is used. If"heterogeneous-circular",
#'   the same as above, but with the constraint that ellipsoids must be spheres.
#'   If "homogeneous-circular", the same as above but the radii of the spheres are
#'   identical. Default is "general". This argument only works for \code{model = "kmeans"}.
#' @param init Methods for choosing initial values of "kmeans" fitting.
#'   Must be "hierarchical" or "kmeans". If "hierarchical", the initial
#'  parameters are obtained with hierarchical clustering method.
#'  If "kmeans", the initial parameters are obtained with extrinsic k-means method.
#'   Additional arguments for k-means clustering and hierarchical clustering can be designated
#'   via argument \code{...}. If no options are designated, \code{nstart=1} for \code{init="kmeans"}
#'   and \code{method="complete"} for \code{init="hierarchical"} are used. Default is "hierarchical".
#' @param d pairwise distance matrix(\code{dist} object) for \code{init = "hierarchical"},
#'   which used in hierarchical clustering. If \code{init = "hierarchical"} and \code{d = NULL},
#'   \code{d} will be automatically filled with \code{ang.pdist(data)}.
#' @param ... Further arguments for argument \code{init}. If \code{init = "kmeans"},
#'   these are for \code{\link[stats]{kmeans}}. If \code{init = "hierarchical"},
#'   there are for \code{\link[stats]{hclust}}.
#' @param additional.condition boolean index.
#'   If \code{TRUE}, a singular matrix will be altered to the scaled identity.
#' @param J A scalar or numeric vector for the number(s) of components for \code{model = c("mixture", "kmeans")}.
#'   Default is \code{J = 4}.
#' @param concentration A scalar or numeric vector for the concentration parameter(s) for \code{model = "kde"}.
#'   Default is \code{concentration = 25}.
#' @param THRESHOLD number for difference between updating and
#'   updated parameters. Default is 1e-10.
#' @param maxiter the maximal number of iteration. Default is 200.
#' @param verbose boolean index, which indicates whether display
#'   additional details as to what the algorithm is doing or
#'   how many loops are done. Moreover, if \code{additional.condition} is
#'   \code{TRUE}, the warning message will be reported.
#' @param kmax the maximal number of kappa. If estimated kappa is
#'   larger than \code{kmax}, then put kappa as \code{kmax}.
#' @return \code{icp.torus} returns an \code{icp.torus} object, containing all values
#'   to compute the conformity score (if \code{J} or \code{concentration} is a
#'   single value). if \code{J} or \code{concentration} is a vector containing
#'   multiple values, then \code{icp.torus} returns a list of \code{icp.torus} objects
#' @export
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
#'
#'   Mardia, K. V., Kent, J. T., Zhang, Z., Taylor, C. C., & Hamelryck, T. (2012). Mixtures of concentrated multivariate sine distributions with applications to bioinformatics. \emph{Journal of Applied Statistics}, 39(11), 2475-2492.
#'
#'   Di Marzio, M., Panzera, A., & Taylor, C. C. (2011). Kernel density estimation on the torus. \emph{Journal of Statistical Planning and Inference}, 141(6), 2156-2173.
#'
#'   Shin, J., Rinaldo, A., & Wasserman, L. (2019). Predictive clustering. \emph{arXiv preprint arXiv:1903.08125}.
#' @examples
#' \donttest{
#' data <- toydata1[, 1:2]
#'
#' icp.torus <- icp.torus(data, model = "kmeans",
#'                        kmeansfitmethod = "general",
#'                        J = 4, concentration = 25)
#' }
icp.torus <- function(data, split.id = NULL,
                      model = c("kmeans", "kde", "mixture"),
                      mixturefitmethod = c("axis-aligned","circular","general"),
                      kmeansfitmethod = c("general", "homogeneous-circular",
                                          "heterogeneous-circular",
                                          "ellipsoids"),
                      init = c("hierarchical", "kmeans"),
                      d = NULL,
                      additional.condition = TRUE,
                      J = 4, concentration = 25, kmax = 500,
                      THRESHOLD = 1e-10, maxiter = 200,
                      verbose = TRUE, ...){
  # returns an icp.torus object, containing all values to compute the conformity score.

  # Use sample splitting to produce (inductive) conformal prediction sets
  # if split.id is not supplied, then use a random set of size floor(n/2) for phat estimation
  # split.id is a n-vector consiting of values 1 (estimation) and 2 (evaluation)

  # param contains the number of components for mixture fitting
  #       and the concentration parameter.

  # if data contains NAs, the rows containing NAs are removed by na.omit()

  data <- stats::na.omit(data)
  if (!is.matrix(data)) {data <- as.matrix(data)}

  # if (is.vector(method)) {method <- "kmeans" }
  # if (is.vector(mixturefitmethod)) {mixturefitmethod <- "axis-aligned" }
  # if (is.vector(kmeansfitmethod)) {kmeansfitmethod <- "general" }
  # if (is.vector(init)){ type <- "hierarchical" }

  data <- on.torus(data)

  model <- match.arg(model)
  mixturefitmethod <- match.arg(mixturefitmethod)
  kmeansfitmethod <- match.arg(kmeansfitmethod)
  init <- match.arg(init)


  if (ncol(data) > 2 && model != "kmeans"){
    warning("kde and mixture methods are not implemented for high dimensional case (>= 3)")
    model <- "kmeans"
  }

  # sample spliting; preparing data
  n <- nrow(data)
  if (is.null(split.id)){
    split.id <- rep(2,n)
    split.id[ sample(n,floor(n/2)) ] <- 1
  }



  # if concentration is a vector, return a list of icp.torus objects
  if(length(concentration) > 1){
    if (model == "kde"){
      icp.torus.objects <- lapply(concentration, function(kappa){
        icp.torus(data, split.id = split.id, model = model,
                        init = init,
                        additional.condition = TRUE,
                        concentration = kappa, ...)
        })
      names(icp.torus.objects) <- paste("conc", concentration, sep = "_")
      return(icp.torus.objects)

    }else{
      concentration = concentration[1]
    }
  }

  # What follows is for single conc. and J.

  X1 <- data[split.id == 1, ]
  X2 <- data[split.id == 2, ]
  n2 <- nrow(X2)

  if (init == "hierarchical" && is.null(d)) {
    d <- ang.pdist(X1)
  }

  # if J is a vector, return a list of icp.torus objects
  if(length(J)>1){
    if(model != "kde"){
      icp.torus.objects <- lapply(J, function(j){
        icp.torus(data, split.id = split.id, model = model,
                  mixturefitmethod = mixturefitmethod,
                  kmeansfitmethod = kmeansfitmethod,
                  init = init,
                  d = d,
                  additional.condition = additional.condition,
                  J = j,
                  kmax = kmax,
                  THRESHOLD = THRESHOLD,
                  maxiter = maxiter,
                  verbose = verbose, ...
                  )
      })
      names(icp.torus.objects) <- paste("J", J, sep = "_")

      return(icp.torus.objects)

    }else{
      J = J[1]
    }
  }

  # Prepare output
  icp.torus <- list(n2 = n2, split.id = split.id, d = ncol(data), data = as.data.frame(data))

  # For each method, use X1 to estimate phat, then use X2 to provide ranks.




  # 1. kde
  if (model == "kde"){
    icp.torus$model <- "kde"

    if (!is.numeric(concentration) | concentration <= 0) {
      concentration <- 25
      warning("Concentration must be a positive number. Reset as concentration = 25 (default)\n")
    }

    phat <- kde.torus(X1, X2, concentration = concentration)
    icp.torus$concentration <- concentration
    icp.torus$score <- sort(phat)
    icp.torus$X1 <- X1

  }

  # 2. mixture fitting by EM
  if (model == "mixture"){
    icp.torus$model <- "mixture"
    icp.torus$fittingmethod <- mixturefitmethod


    if (!is.numeric(J) | J < 1) {
      J <- 4
      warning("The number of components must be a positive integer Reset as J = 4 (default)\n")
    }
    J = as.integer(J)

    vm2mixfit <- EMsinvMmix(X1, J = J, parammat = EMsinvMmix.init(data, J),
                            THRESHOLD = THRESHOLD, maxiter = maxiter,
                            type = mixturefitmethod,
                            kmax = kmax,
                            verbose = verbose)
    icp.torus$fit <- vm2mixfit

    # compute phat(X2)
    phat <- BAMBI::dvmsinmix(X2,kappa1 = vm2mixfit$parammat[2, ],
                             kappa2 = vm2mixfit$parammat[3, ],
                             kappa3 = vm2mixfit$parammat[4, ],
                             mu1 = vm2mixfit$parammat[5, ],
                             mu2 = vm2mixfit$parammat[6, ],
                             pmix = vm2mixfit$parammat[1, ], log = FALSE)
    icp.torus$score <- sort(phat)


    # compute phat_max(X2)
    phatj <- phat.eval(X2, vm2mixfit$parammat)
    icp.torus$score_max <- sort(do.call(pmax, as.data.frame(phatj)))

    # compute parameters for ellipses and phat_e(X2)
    ellipse.param <- norm.appr.param(vm2mixfit$parammat)
    icp.torus$ellipsefit <- ellipse.param

    ehatj <- ehat.eval(X2, ellipse.param)
    icp.torus$score_ellipse <- sort(do.call(pmax, as.data.frame(ehatj)))

  }

  # 3. elliptical kmeans algorithm
  if (model == "kmeans"){
    # implement extrinsic kmeans clustering for find the centers
    icp.torus$model <- "kmeans"
    icp.torus$fittingmethod <- kmeansfitmethod

    if (!is.numeric(J) | J < 1) {
      J <- 4
      warning("The number of components must be a positive integer Reset as J = 4 (default)\n")
    }
    J = as.integer(J)

    # consider -R as ehat in von mises mixture approximation
    # where R is the notation in J. Shin (2019)
    ellipse.param <- ellip.kmeans.torus(X1, centers = J,
                                        type = kmeansfitmethod,
                                        init = init,
                                        additional.condition = additional.condition,
                                        THRESHOLD = THRESHOLD, maxiter = maxiter,
                                        verbose = verbose, d = d, ...)

    icp.torus$ellipsefit <- ellipse.param

    ellipsej <- ehat.eval(X2, ellipse.param)
    icp.torus$score_ellipse <- sort(do.call(pmax, as.data.frame(ellipsej)))

  }
  return(structure(icp.torus, class = "icp.torus"))

}

#' Inductive prediction sets for each level
#'
#' \code{icp.torus.eval} evaluates whether each pre-specified evaluation point
#'   is contained in the inductive conformal prediction sets for each given
#'   level.
#'
#' @param icp.torus an object containing all values to compute the conformity
#'   score, which will be constructed with \code{icp.torus}.
#' @param level either a scalar or a vector, or even \code{NULL}. Default value
#'   is 0.1.
#' @param eval.point N x N numeric matrix on \eqn{[0, 2\pi)^2}. Default input is
#'   \code{grid.torus}.
#' @return returns a \code{cp} object with the boolean values which
#'   indicate whether each evaluation point is contained in the inductive
#'   conformal prediction sets for each given level.
#' @export
#' @seealso \code{\link[ClusTorus]{grid.torus}}, \code{\link[ClusTorus]{icp.torus}}
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
#' @examples
#' \donttest{
#' data <- toydata1[, 1:2]
#'
#' icp.torus <- icp.torus(data, model = "kmeans",
#'                        mixturefitmethod = "general",
#'                        J = 4, concentration = 25)
#'
#' icp.torus.eval(icp.torus, level = c(0.1, 0.08), eval.point = grid.torus())
#' }
icp.torus.eval <- function(icp.torus, level = 0.1, eval.point = grid.torus()){
  # evaluates Chat_kde, Chat_mix, Chat_max, Chat_e.
  stopifnot(class(icp.torus) == "icp.torus")
  N <- nrow(eval.point)

  n2 <- icp.torus$n2
  nalpha <- length(level)
  cp <- list(level = level, eval.point = eval.point)

  ialpha <- floor((n2 + 1) * level)

  if(icp.torus$model == "kde"){
    cp$model <- "kde"
    Chat_kde <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_kde) <- level

    phat.grid <- kde.torus(icp.torus$X1, eval.point, concentration = icp.torus$concentration)

    # for (i in 1:nalpha){
    #   ialpha <- floor((n2 + 1) * level[i])
    #
    #   # indices for inclusion in Chat_kde
    #   Chat_kde[, i] <- phat.grid >= icp.torus$kde$score[ialpha]
    # }
    scores_kde <- icp.torus$score[ialpha]
    Chat_kde <- sweep(replicate(nalpha, phat.grid), 2, scores_kde, ">=")

    cp$Chat_kde <- Chat_kde
  }

  if(icp.torus$model == "mixture"){
    cp$model <- "mixture"
    Chat_mix <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_mix) <- level

    Chat_max <- Chat_mix
    Chat_e <- Chat_mix

    phatj <- phat.eval(eval.point, icp.torus$fit$parammat)
    ehatj <- ehat.eval(eval.point, icp.torus$ellipsefit)
    phat_mix <- rowSums(phatj)
    # phat_max <- apply(phatj, 1, max)
    # ehat <- apply(ehatj, 1, max)
    phat_max <- do.call(pmax, as.data.frame(phatj))
    ehat <- do.call(pmax, as.data.frame(ehatj))


    # for (i in 1:nalpha){
    #   ialpha <- floor((n2 + 1) * level[i])
    #
    #   Chat_mix[, i] <- phat_mix >= icp.torus$mixture$score[ialpha]
    #   Chat_max[, i] <- phat_max >= icp.torus$mixture$score_max[ialpha]
    #   Chat_e[, i]   <-    ehat  >= icp.torus$mixture$score_ellipse[ialpha]
    # }
    scores_mix <- icp.torus$score[ialpha]
    scores_max <- icp.torus$score_max[ialpha]
    scores_e <- icp.torus$score_ellipse[ialpha]

    Chat_mix <- sweep(replicate(nalpha, phat_mix), 2, scores_mix, ">=")
    Chat_max <- sweep(replicate(nalpha, phat_max), 2, scores_max, ">=")
    Chat_e <- sweep(replicate(nalpha, ehat), 2, scores_e, ">=")

    cp$Chat_mix <- Chat_mix
    cp$Chat_max <- Chat_max
    cp$Chat_e <- Chat_e
  }

  if(icp.torus$model == "kmeans"){
    cp$model <- "kmeans"
    Chat_kmeans <- matrix(0, nrow = N, ncol = nalpha)

    ellipsej <- ehat.eval(eval.point, icp.torus$ellipsefit)
    ellipse <- do.call(pmax, as.data.frame(ellipsej))

    # for (i in 1:nalpha){
    #   ialpha <- floor((n2 + 1) * level[i])
    #   Chat_kmeans[, i] <- sphere >= icp.torus$kmeans$score_sphere[ialpha]
    # }
    scores_ellipse <- icp.torus$score_ellipse[ialpha]
    Chat_kmeans <- sweep(replicate(nalpha, ellipse), 2, scores_ellipse, ">=")

    cp$Chat_kmeans <- Chat_kmeans

  }

  return(structure(cp, class = "icp.torus.eval"))
}
