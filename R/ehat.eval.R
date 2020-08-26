#' Normal approximation to the log - weighted bivariate von Mises sine density
#'
#' \code{ehat.eval} evaluates the approximation of
#'   log - weighted bivariate von Mises for each given data and f
#'   or each given parameters.
#'
#' @param X n x 2 toroidal data on \eqn{[-\pi, \pi)^2}
#' @param ellipse.param list which is consisting of mean of each angular
#'   coordinate, inverse of each covariance matrix, and constant term
#' @return nrow(X) times ncol(parammat) (n x J) matrix
#' @export

ehat.eval <- function(X, ellipse.param){

  # evaluate ehat_j(x). Returns nrow(X) times ncol(parammat) (n x J) matrix.
  # ehat(x) is the columnwise maximum of ehat_j(x): apply(ehatj,1,max)
  # used for elliptically-approximated mixture models
  n2 <- nrow(X)
  J <- length(ellipse.param$mu1)

  ehatj <- matrix(0, nrow = n2, ncol = J)
  for(j in 1:J){
    z <- tor.minus(X, c(ellipse.param$mu1[j], ellipse.param$mu2[j]))
    S <- ellipse.param$Sigmainv[[j]]
    A <- z %*% S
    ehatj[,j] <- -apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]}) + ellipse.param$c[j]
  }
  ehatj
}
