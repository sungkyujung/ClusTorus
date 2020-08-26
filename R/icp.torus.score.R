#' Inductive conformal prediction sets
#'
#' \code{icp.torus.score} returns an icp.torus object, containing all values
#'   to compute the conformity score.
#'
#' @param data n x 2 matrix of toroidal data on \eqn{[-\pi, \pi)^2}
#' @param split.id a n-vector consisting of values 1 (estimation) and 2(evaluation)
#' @param method a string one of "all", "kde" and "mixture" which determines the
#'   phat.
#' @param mixturefitmethod a string one of "circular", "axis-aligned", "general",
#'   and "Bayesian" which determines the fitting method.("Bayesian" is not yet
#'   supported)
#' @param param the number of components (in \code{list} form) for mixture
#'   fitting and the concetnration parameter.
#' @return returns an icp.torus object, containing all values
#'   to compute the conformity score
#' @export
#' @seealso \code{\link{EMsinvMmix}}, \code{\link[BAMBI]{dvmsinmix}}
#' @examples
#' \dontrun{
#' data <- matrix(c(-pi/3, -pi/3, pi/2, pi/4),
#'                nrow = 2, byrow = TRUE)
#'
#' icp.torus.score(data, method = "kde", mixturefitmethod = "circular")
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
