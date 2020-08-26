#' Mahalanobis distance to ellipse
#'
#' \code{mah.dist2.eval} evaluates the mahalanobis distance from a toroidal
#'   point to each ellipses.
#'
#' @inheritParams ehat.eval
#' @return \code{mah.dist2.eval} returns  nrow(X) times ncol(parammat)
#'   (n x J) matrix.
#' @export
mah.dist2.eval <- function(X, ellipse.param){

  # evaluate mahalanobis_distance from x to each ellipses. Returns nrow(X) times ncol(parammat) (n x J) matrix.

  n2 <- nrow(X)
  J <- length(ellipse.param$mu1)

  ehatj <- matrix(0,nrow = n2,ncol = J)
  for(j in 1:J){
    z <- tor.minus(X, c(ellipse.param$mu1[j], ellipse.param$mu2[j]) )
    S <- ellipse.param$Sigmainv[[j]]
    A <- z %*% S
    ehatj[,j] <-  apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]})
  }
  ehatj
}
