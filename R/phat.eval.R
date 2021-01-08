# Product of weights and von Mises sine distribution functions
#
# \code{phat.eval} evaluates \eqn{pi_j * f_j(x)}, where \eqn{pi_j} is a weight
#   and \eqn{f_j(x)} is the j-th von Mises sine density which is determined by
#   j-th parameters of the given parameter data.
#
# @param X n x 2 toroidal data on \eqn{[0, 2\pi)^2}
# @param parammat 6 x J parameter data with the following components:
#
#   \code{parammat[1, ]} : the weights for each von Mises sine density
#
#   \code{parammat[n + 1, ]} : \eqn{\kappa_n} for each von Mises
#   sine density for n = 1, 2, 3
#
#   \code{parammat[m + 4 , ]} : \eqn{\mu_m} for each von Mises
#   sine density for m = 1, 2
# @return returns nrow(X) times ncol(parammat) (n x J) matrix
# @seealso \code{\link[BAMBI]{dvmsin}}
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
#                      0, 2, 4), nrow = 6, byrow = TRUE)
#
# phat.eval(X, parammat)
# }

phat.eval <- function(X, parammat){
  # evaluate pi_j f_j(x). Returns nrow(X) times ncol(parammat) (n x J) matrix.
  # used for mixture models
  n2 <- nrow(X)
  J <- ncol(parammat)
  phatj <- matrix(0,nrow = n2,ncol = J)
  for(j in 1:J){
    phatj[, j] <- BAMBI::dvmsin(X, kappa1 = parammat[2, j],
                               kappa2 = parammat[3, j],
                               kappa3 = parammat[4, j],
                               mu1 = parammat[5, j],
                               mu2 = parammat[6, j], log = FALSE
    ) * parammat[1, j]
  }
  phatj <- ifelse(is.nan(phatj), 0, phatj)
  return(phatj)
}
