#' Conformity score for inductive prediction sets
#'
#' \code{icp.torus.score} prepares all values
#'   for computing the conformity score for specified methods.
#'
#' @param data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#'   or \eqn{[-\pi, \pi)^d}
#' @param split.id a n-dimensinal vector consisting of values 1 (estimation)
#'   and 2(evaluation)
#' @param method a string one of "all", "kde", "mixture", and "kmeans"
#'   which determines the model for clustering. Default is "all". Moreover, if
#'   the dimension of data space is larger than 2, then automatically
#'   "kmeans".
#' @param mixturefitmethod a string one of "circular", "axis-aligned", "general",
#'   and "Bayesian" which determines the fitting mixture method.("Bayesian" is not yet
#'   supported) Default is "axis-aligned".
#' @param kmeansfitmethod character which must be "homogeneous-circular",
#'  "heterogeneous-circular", or "general".
#'   If "homogeneous-circular", the radii of k-spheres are identical.
#'   If "heterogeneous-cricular", the radii of k-spheres may be different.
#'   If "ellipsoids", cluster with k-ellipsoids without optimized parameters.
#'   If, "general", clustering with k-ellipses. The parameters to construct
#'   the ellipses are optimized with generalized Lloyd algorithm, which is
#'   modified for toroidal space. To see the detail, see the references.
#' @param init determine the initial parameter of "kmeans" method,
#'   for option "general". Must be "kmeans" or "hierarchical".
#'   If "kmeans", the initial parameters are obtained with extrinsic kmeans
#'   method.
#'   If "hierarchical", the initial parameters are obtained with hierarchical
#'   clustering method.
#'   Default is "kmeans".
#' @param additional.condition boolean index.
#'   If \code{TRUE}, a singular matrix will be altered to the scaled identity.
#' @param param the number of components (in \code{list} form) for mixture
#'   fitting and the concetnration parameter.
#' @param THRESHOLD number for difference between updating and
#'   updated parameters. Default is 1e-10.
#' @param maxiter the maximal number of iteration.
#' @param verbose boolean index, which indicates whether display
#'   additional details as to what the algorithm is doing or
#'   how many loops are done. Moreover, if \code{additional.condition} is
#'   \code{TRUE}, the warning message will be reported.
#' @param kmax the maximal number of kappa. If estimated kappa is
#'   larger than \code{kmax}, then put kappa as \code{kmax}.
#' @return returns an \code{icp.torus} object, containing all values
#'   to compute the conformity score.
#' @export
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#' @examples
#' \donttest{
#' data <- ILE[1:200, 1:2]
#'
#' icp.torus <- icp.torus.score(data, method = "all",
#'                              mixturefitmethod = "general",
#'                              kmeansfitmethod = "general",
#'                              param = list(J = 4, concentration = 25))
#' }
icp.torus.score <- function(data, split.id = NULL,
                            method = c("all", "kde", "mixture", "kmeans"),
                            mixturefitmethod = c("circular", "axis-aligned", "general", "Bayesian"),
                            kmeansfitmethod = c("homogeneous-circular",
                                                "heterogeneous-circular",
                                                "ellipsoids",
                                                "general"),
                            init = c("kmeans", "hierarchical"),
                            additional.condition = TRUE,
                            param = list(J = 4, concentration = 25), kmax = 500,
                            THRESHOLD = 1e-10, maxiter = 200,
                            verbose = TRUE){
  # returns an icp.torus object, containing all values to compute the conformity score.

  # Use sample splitting to produce (inductive) conformal prediction sets
  # if split.id is not supplied, then use a random set of size floor(n/2) for phat estimation
  # split.id is a n-vector consiting of values 1 (estimation) and 2 (evaluation)

  # param contains the number of components for mixture fitting
  #       and the concentration parameter.
  if (!is.matrix(data)) {data <- as.matrix(data)}

  if (is.null(method)) {method <- "all" }
  if (is.null(mixturefitmethod)) {mixturefitmethod <- "axis-aligned" }
  if (is.null(kmeansfitmethod)) {kmeansfitmethod <- "omogeneous-circular" }
  if (is.null(init)){ type <- "kmeans" }

  data <- on.torus(data)

  method <- match.arg(method)
  mixfitmethod <- match.arg(mixturefitmethod)
  kmeansfitmethod <- match.arg(kmeansfitmethod)

  if (ncol(data) > 2 && method != "kmeans"){
    warning("kde and mixture methods are not implemented for high dimensional case (>= 3)")
    method <- "kmeans"
  }

  # sample spliting; preparing data
  n <- nrow(data)
  if (is.null(split.id)){
    split.id <- rep(2,n)
    split.id[ sample(n,floor(n/2)) ] <- 1
  }

  X1 <- data[split.id == 1, ]
  X2 <- data[split.id == 2, ]
  n2 <- nrow(X2)

  # Prepare output
  icp.torus <- list(kde = NULL, mixture = NULL, kmeans = NULL,
                    n2 = n2, split.id = split.id, d = ncol(data))

  # For each method, use X1 to estimate phat, then use X2 to provide ranks.

  # 1. kde
  if (sum(method == c("kde", "all")) == 1){

    phat <- kde.torus(X1, X2, concentration = param$concentration)
    # phat.X2.sorted <- sort(phat)

    icp.torus$kde$concentration <- param$concentration
    icp.torus$kde$score <- sort(phat)
    icp.torus$kde$X1 <- X1

  }

  # 2. mixture fitting
  if (sum(method == c("mixture", "all")) == 1){

    icp.torus$mixture$fittingmethod <- mixturefitmethod

    if (mixturefitmethod != "Bayesian"){
      vm2mixfit <- EMsinvMmix(X1, J = param$J, parammat = EMsinvMmix.init(data, param$J),
                              THRESHOLD = THRESHOLD, maxiter = maxiter,
                              type = mixturefitmethod,
                              kmax = kmax,
                              verbose = verbose)
      icp.torus$mixture$fit <- vm2mixfit
    } else {
      #vm2mixfit ## USE BAMBI
      stop("Bayesian not yet implemented")
    }

    # compute phat(X2)
    phat <- BAMBI::dvmsinmix(X2,kappa1 = vm2mixfit$parammat[2, ],
                             kappa2 = vm2mixfit$parammat[3, ],
                             kappa3 = vm2mixfit$parammat[4, ],
                             mu1 = vm2mixfit$parammat[5, ],
                             mu2 = vm2mixfit$parammat[6, ],
                             pmix = vm2mixfit$parammat[1, ], log = FALSE)
    icp.torus$mixture$score <- sort(phat)


    # compute phat_max(X2)
    phatj <- phat.eval(X2, vm2mixfit$parammat)

    # phatj <- matrix(0,nrow = n2,ncol = param$J)
    # for(j in 1:param$J){
    #   phatj[,j] <- BAMBI::dvmsin(X2, kappa1 = vm2mixfit$parammat[2,j],
    #                              kappa2 = vm2mixfit$parammat[3,j],
    #                              kappa3 = vm2mixfit$parammat[4,j],
    #                              mu1 = vm2mixfit$parammat[5,j],
    #                              mu2 = vm2mixfit$parammat[6,j],log = FALSE
    #   ) * vm2mixfit$parammat[1,j]
    # }
    # icp.torus$mixture$score_max <- sort(apply(phatj, 1, max))
    icp.torus$mixture$score_max <- sort(do.call(pmax, as.data.frame(phatj)))

    # compute parameters for ellipses and phat_e(X2)
    ellipse.param <- norm.appr.param(vm2mixfit$parammat)
    # ellipse.param <- list(mu1 = NULL, mu2=NULL, Sigmainv = NULL,c = NULL)
    # ellipse.param$mu1 <- vm2mixfit$parammat[5,]
    # ellipse.param$mu2 <- vm2mixfit$parammat[6,]
    # for (j in 1:param$J){
    #   kap1 <- vm2mixfit$parammat[2,j]
    #   kap2 <- vm2mixfit$parammat[3,j]
    #   lamb <- vm2mixfit$parammat[4,j]
    #   pi_j <- vm2mixfit$parammat[1,j]
    #   ellipse.param$Sigmainv[[j]] <-
    #     matrix(c(kap1, rep(lamb,2), kap2), nrow = 2)
    #   ellipse.param$c[j] <- 2*log(pi_j * (kap1*kap2 - lamb^2) )
    # }
    icp.torus$mixture$ellipsefit <- ellipse.param

    ehatj <- ehat.eval(X2, ellipse.param)
    # ehatj <- matrix(0,nrow = n2,ncol = param$J)
    # for(j in 1:param$J){
    #   z <- tor.minus(X2, c(ellipse.param$mu1[j], ellipse.param$mu2[j]) )
    #   S <- ellipse.param$Sigmainv[[j]]
    #   A <- z %*% S
    #   ehatj[,j] <- -apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]}) + ellipse.param$c[j]
    # }
    # icp.torus$mixture$score_ellipse <- sort(apply(ehatj, 1, max))
    icp.torus$mixture$score_ellipse <- sort(do.call(pmax, as.data.frame(ehatj)))

  }

  # 3. kmeans to kspheres
  if (sum(method == c("kmeans", "all")) == 1){
    # implement extrinsic kmeans clustering for find the centers

    # consider -R as ehat in von mises mixture approximation
    # where R is the notation in J. Shin (2019)
    sphere.param <- kmeans.kspheres(X1, centers = param$J,
                                    type = kmeansfitmethod,
                                    init = init,
                                    additional.condition = additional.condition,
                                    THRESHOLD = THRESHOLD, maxiter = maxiter,
                                    verbose = verbose)

    icp.torus$kmeans$spherefit <- sphere.param

    spherej <- ehat.eval(X2, sphere.param)
    # icp.torus$kmeans$score_sphere <- sort(apply(spherej, 1, max))
    icp.torus$kmeans$score_sphere <- sort(do.call(pmax, as.data.frame(spherej)))

  }

  return(icp.torus)

}

#' Inductive prediction sets for each level
#'
#' \code{icp.torus.eval} evaluates whether each pre-specified evaluation point
#'   is contained in the inductive conformal prediction sets for each given
#'   level.
#'
#' @param icp.torus an object containing all values to compute the conformity
#'   score, which will be constructed with \code{icp.torus.score}.
#' @param level either a scalar or a vector, or even \code{NULL}. Default value
#'   is 0.1.
#' @param eval.point N x N numeric matrix on \eqn{[0, 2\pi)^2}. Default input is
#'   \code{grid.torus}.
#' @return returns a \code{cp} object with the boolean values which
#'   indicate whether each evaluation point is contained in the inductive
#'   conformal prediction sets for each given level.
#' @export
#' @seealso \code{\link[ClusTorus]{grid.torus}}, \code{\link[ClusTorus]{icp.torus.score}}
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#' @examples
#' \donttest{
#' data <- ILE[1:200, 1:2]
#'
#' icp.torus <- icp.torus.score(data, method = "all",
#'                              mixturefitmethod = "general",
#'                              param = list(J = 4, concentration = 25))
#'
#' icp.torus.eval(icp.torus, level = c(0.1, 0.08), eval.point = grid.torus())
#' }

icp.torus.eval <- function(icp.torus, level = 0.1, eval.point = grid.torus()){
  # evaluates Chat_kde, Chat_mix, Chat_max, Chat_e.
  N <- nrow(eval.point)

  n2 <- icp.torus$n2
  nalpha <- length(level)
  cp <- list(Chat_kde = NULL, Chat_mix = NULL, Chat_max = NULL, Chat_e = NULL,
             Chat_kmeans = NULL,
             level = level,
             eval.point = eval.point)

  ialpha <- floor((n2 + 1) * level)


  if(!is.null(icp.torus$kde)){

    Chat_kde <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_kde) <- level

    phat.grid <- kde.torus(icp.torus$kde$X1, eval.point, concentration = icp.torus$kde$concentration)

    # for (i in 1:nalpha){
    #   ialpha <- floor((n2 + 1) * level[i])
    #
    #   # indices for inclusion in Chat_kde
    #   Chat_kde[, i] <- phat.grid >= icp.torus$kde$score[ialpha]
    # }
    scores_kde <- icp.torus$kde$score[ialpha]
    Chat_kde <- sweep(replicate(nalpha, phat.grid), 2, scores_kde, ">=")

    cp$Chat_kde <- Chat_kde
  }

  if(!is.null(icp.torus$mixture)){
    Chat_mix <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_mix) <- level

    Chat_max <- Chat_mix
    Chat_e <- Chat_mix

    phatj <- phat.eval(eval.point, icp.torus$mixture$fit$parammat)
    ehatj <- ehat.eval(eval.point, icp.torus$mixture$ellipsefit)
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
    scores_mix <- icp.torus$mixture$score[ialpha]
    scores_max <- icp.torus$mixture$score_max[ialpha]
    scores_e <- icp.torus$mixture$score_ellipse[ialpha]

    Chat_mix <- sweep(replicate(nalpha, phat_mix), 2, scores_mix, ">=")
    Chat_max <- sweep(replicate(nalpha, phat_max), 2, scores_max, ">=")
    Chat_e <- sweep(replicate(nalpha, ehat), 2, scores_e, ">=")

    cp$Chat_mix <- Chat_mix
    cp$Chat_max <- Chat_max
    cp$Chat_e <- Chat_e
  }

  if(!is.null(icp.torus$kmeans)){
    Chat_kmeans <- matrix(0, nrow = N, ncol = nalpha)

    spherej <- ehat.eval(eval.point, icp.torus$kmeans$spherefit)
    sphere <- do.call(pmax, as.data.frame(spherej))

    # for (i in 1:nalpha){
    #   ialpha <- floor((n2 + 1) * level[i])
    #   Chat_kmeans[, i] <- sphere >= icp.torus$kmeans$score_sphere[ialpha]
    # }
    scores_sphere <- icp.torus$kmeans$score_sphere[ialpha]
    Chat_kmeans <- sweep(replicate(nalpha, sphere), 2, scores_sphere, ">=")

    cp$Chat_kmeans <- Chat_kmeans

  }

  return(cp)


}
