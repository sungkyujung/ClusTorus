# Elliptical approximation of bivariate von Mises
#
# \code{norm.appr.param} approximates the bivariate von Mises sine densities
#   to bivariate normal densities with returning the parameters for
#   the normal densities.
# @param parammat 6 x J parameter data with the following components:
#
#   \code{parammat[1, ]} : the weights for each von Mises sine density
#
#   \code{parammat[n + 1, ]} : \eqn{\kappa_n} for each von Mises
#   sine density for n = 1, 2, 3
#
#   \code{parammat[m + 4, ]} : \eqn{\mu_m} for each von Mises
#    sine density for m = 1, 2
# @return returns approximated parameters for bivariate normal
# distribution with \code{list}:
#
#   \code{list$Sigmainv[j]} : approximated covariance matrix for
#   j-th bivariate normal distribution, approximation of the j-th von Mises.
#
#   \code{list$c[j]} : approximated \eqn{|2\pi\Sigma|^-1} for
#   j-th bivariate normal distribution, approximation of the j-th von Mises.
# @seealso \code{\link[BAMBI]{dvmsin}}
# @examples
# \dontrun{
#
# parammat <- matrix(c(0.4, 0.3, 0.3,
#                      20, 25, 25,
#                      30, 25, 20,
#                      1, 2, 3,
#                      1, 2, 3,
#                      0, 2, 4), nrow = 6, byrow =TRUE)
#
# norm.appr.param(parammat)
# }
norm.appr.param <- function(parammat){
  # parameters for elliptical approximation of bivariate von mises
  J <- ncol(parammat)

  ellipse.param <- list(mu = NULL, Sigmainv = NULL, c = NULL)
  ellipse.param$mu <- cbind(parammat[5, ], parammat[6, ])
  for (j in 1:J){
    kap1 <- parammat[2, j]
    kap2 <- parammat[3, j]
    lamb <- parammat[4, j]
    pi_j <- parammat[1, j]

    # # small adjustment for low concentration
    # if (kap1 < 10){
    #   lamb <- lamb * (kap1/10)^(1/10)
    #   kap1 <- kap1 * (kap1/10)^(1/10)
    # }
    # if (kap2 < 10){
    #   lamb <- lamb * (kap2/10)^(1/10)
    #   kap2 <- kap1 * (kap2/10)^(1/10)
    # }
    #

    ellipse.param$Sigmainv[[j]] <-
      matrix(c(kap1, rep(lamb, 2), kap2), nrow = 2)
    ellipse.param$c[j] <- log(pi_j^2 * (kap1 * kap2 - lamb^2))
  }

  return(ellipse.param)
}
