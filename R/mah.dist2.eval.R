# Mahalanobis distance to ellipse
#
# \code{mah.dist2.eval} evaluates the mahalanobis distance from a toroidal
#   point to each ellipses.
#
# @inheritParams ehat.eval
# @return \code{mah.dist2.eval} returns nrow(X) times ncol(parammat)
#   (n x J) matrix whose entries implies the Mahalanobis distance
#   between i-th data point and the j-th ellipse.
# @examples
# \dontrun{
# ## mean vectors
#
# Mu1 <- c(3, 0)
# Mu2 <- c(2, 2)
# Mu3 <- c(1, 4)
#
# ## covariance matrices
#
# Sigma1 <- matrix(c(0.1, 0.05, 0.05, 0.2), 2, 2)
# Sigma2 <- matrix(c(0.1, 0, 0, 0.01), 2, 2)
# Sigma3 <- matrix(c(0.01, 0, 0, 0.1), 2, 2)
#
# ## 2-dimensional multivariate normal data wrapped with toroidal space
# require(MASS)
# X <- rbind(mvrnorm(n=70, Mu1, Sigma1),
#            mvrnorm(n=50, Mu2, Sigma2),
#            mvrnorm(n=50, Mu3, Sigma3))
#
# X <- on.torus(X)
#
# parammat <- matrix(c(0.4, 0.3, 0.3,
#                      20, 25, 25,
#                      30, 25, 20,
#                      1, 2, 3,
#                      1, 2, 3,
#                      0, 2, 4), nrow = 6, byrow =TRUE)
#
# elipse.param <- norm.appr.param(parammat)
#
# mah.dist2.eval(X, ellipse.param)
# }
mah.dist2.eval <- function(X, ellipse.param){

  # evaluate mahalanobis_distance from x to each ellipses. Returns nrow(X) times ncol(parammat) (n x J) matrix.

  n2 <- nrow(X)
  J <- length(ellipse.param$c)

  ehatj <- matrix(0,nrow = n2,ncol = J)
  for(j in 1:J){
    # z <- tor.minus(X, c(ellipse.param$mu1[j], ellipse.param$mu2[j]) )
    z <- tor.minus(X, ellipse.param$mu[j, ])
    S <- ellipse.param$Sigmainv[[j]]
    A <- z %*% S
    # ehatj[,j] <-  apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]})

    ehatj[,j] <- rowSums(A * z)
  }
  ehatj
}
