# Normal approximation to the log - weighted bivariate von Mises sine density
#
# \code{ehat.eval} evaluates the approximation of
#   log - weighted bivariate von Mises for each given data and f
#   or each given parameters.
#
# @param X n x 2 toroidal data on \eqn{[0, 2\pi)^2}
# @param ellipse.param list which is consisting of mean of each angular
#   coordinate, inverse of each covariance matrix, and constant term
# @return nrow(X) times ncol(parammat) (n x J) matrix
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
# ellipse.param <- norm.appr.param(parammat)
#
# ehat.eval(X, ellipse.param)
# }

ehat.eval <- function(X, ellipse.param){

  # evaluate ehat_j(x). Returns nrow(X) times ncol(parammat) (n x J) matrix.
  # ehat(x) is the columnwise maximum of ehat_j(x): apply(ehatj,1,max)
  # used for elliptically-approximated mixture models
  n2 <- nrow(X)
  J <- length(ellipse.param$c)

  ehatj <- matrix(0, nrow = n2, ncol = J)
  for(j in 1:J){
    # z <- tor.minus(X, c(ellipse.param$mu1[j], ellipse.param$mu2[j]))
    z <- tor.minus(X, ellipse.param$mu[j, ])
    A <- z %*% ellipse.param$Sigmainv[[j]]
    # ehatj[,j] <- -apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]}) + ellipse.param$c[j]
    # ehatj[,j] <- -diag(A %*% t(z)) + ellipse.param$c[j]
    ehatj[,j] <- -rowSums(A * z) + ellipse.param$c[j]
  }
  ehatj
}
