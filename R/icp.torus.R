#' Conformity score for inductive prediction sets
#'
#' \code{icp.torus.score} returns an icp.torus object, containing all values
#'   to compute the conformity score.
#'
#' @param data n x 2 matrix of toroidal data on \eqn{[0, 2\pi)^2}
#' @param split.id a n-dimensinal vector consisting of values 1 (estimation)
#'   and 2(evaluation)
#' @param method a string one of "all", "kde" and "mixture" which determines the
#'   phat.
#' @param mixturefitmethod a string one of "circular", "axis-aligned", "general",
#'   and "Bayesian" which determines the fitting method.("Bayesian" is not yet
#'   supported)
#' @param param the number of components (in \code{list} form) for mixture
#'   fitting and the concetnration parameter.
#' @return returns an \code{icp.torus} object, containing all values
#'   to compute the conformity score.
#' @export
#' @seealso \code{\link{EMsinvMmix}}, \code{\link[BAMBI]{dvmsinmix}}
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
#'                              mixturefitmethod = "general",)
#'                              param = list(J = 4, concentration = 25))
#' }
icp.torus.score <- function(data, split.id = NULL,
                            method = c("all", "kde", "mixture"),
                            mixturefitmethod = c("circular", "axis-aligned", "general", "Bayesian"),
                            param = list(J = 4, concentration = 25)){
  # returns an icp.torus object, containing all values to compute the conformity score.

  # Use sample splitting to produce (inductive) conformal prediction sets
  # if split.id is not supplied, then use a random set of size floor(n/2) for phat estimation
  # split.id is a n-vector consiting of values 1 (estimation) and 2 (evaluation)

  # param contains the number of components for mixture fitting
  #       and the concentration parameter.

  method <- match.arg(method)
  mixfitmethod <- match.arg(mixturefitmethod)


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
  icp.torus <- list(kde = NULL, mixture = NULL, n2 = n2, split.id = split.id)

  # For each method, use X1 to estimate phat, then use X2 to provide ranks.


  # 1. kde
  if (method != "mixture"){
    phat <- kde.torus(X1, X2, concentration = param$concentration)
    # phat.X2.sorted <- sort(phat)

    icp.torus$kde$concentration <- param$concentration
    icp.torus$kde$score <- sort(phat)
    icp.torus$kde$X1 <- X1
  }

  # 2. mixture fitting
  if (method != "kde"){

    icp.torus$mixture$fittingmethod <- mixturefitmethod

    if (mixturefitmethod != "Bayesian"){
      vm2mixfit <- EMsinvMmix(X1, J = param$J, parammat = EMsinvMmix.init(data, param$J),
                              THRESHOLD = 1e-10,
                              type = mixturefitmethod,
                              kmax = 500,
                              verbose = TRUE)
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
    icp.torus$mixture$score_max <- sort(apply(phatj, 1, max))

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
    icp.torus$mixture$score_ellipse <- sort(apply(ehatj, 1, max))

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
#' @seealso \code{\link{grid.torus}}, \code{\link{icp.torus.score}}
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
#'                              mixturefitmethod = "general",)
#'                              param = list(J = 4, concentration = 25)
#'
#' icp.torus.eval(icp.torus, level = c(0.1, 0.08), eval.point = grid.torus())
#' }

icp.torus.eval <- function(icp.torus, level = 0.1, eval.point = grid.torus()){
  # evaluates Chat_kde, Chat_mix, Chat_max, Chat_e.
  N <- nrow(eval.point)

  n2 <- icp.torus$n2
  nalpha <- length(level)
  cp <- list(Chat_kde = NULL, Chat_mix = NULL, Chat_max = NULL, Chat_e = NULL,
             level = level,
             phi = eval.point[, 1],
             psi = eval.point[, 2])



  if(!is.null(icp.torus$kde)){

    Chat_kde <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_kde) <- level

    phat.grid <- kde.torus(icp.torus$kde$X1, eval.point, concentration = icp.torus$kde$concentration)
    for (i in 1:nalpha){
      ialpha <- floor((n2 + 1) * level[i])

      # indices for inclusion in Chat_kde
      Chat_kde[, i] <- phat.grid >= icp.torus$kde$score[ialpha]
    }
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
    phat_max <- apply(phatj, 1, max)
    ehat <- apply(ehatj, 1, max)

    for (i in 1:nalpha){
      ialpha <- floor( (n2 + 1) * level[i])

      Chat_mix[, i] <- phat_mix >= icp.torus$mixture$score[ialpha]
      Chat_max[, i] <- phat_max >= icp.torus$mixture$score_max[ialpha]
      Chat_e[, i]   <-    ehat  >= icp.torus$mixture$score_ellipse[ialpha]
    }

    cp$Chat_mix <- Chat_mix
    cp$Chat_max <- Chat_max
    cp$Chat_e <- Chat_e
  }

  return(cp)

}
