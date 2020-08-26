#' Product of weights and von Mises sine distribution functions
#'
#' \code{phat.eval} evaluates \eqn{pi_j * f_j(x)}, where \eqn{pi_j} is a weight
#'   and \eqn{f_j(x)} is the j-th von Mises sine density which is determined by
#'   j-th parameters of the given parameter data.
#'
#' @param X n x 2 toroidal data on \eqn{[-\pi, pi)^2}
#' @param parammat 6 x J parameter data with the following components:
#'   \code{parammat[1, ]} : the weights for each von Mises sine density
#'   \code{parammat[n + 1, ]} : \eqn{\kappa_n} for each von Mises sine density for n = 1, 2, 3
#'   \code{parammat[m + 4 , ]} : \eqn{\mu_m} for each von Mises sine density for m = 1, 2
#' @return returns nrow(X) times ncol(parammat) (n x J) matrix
#' @export
#' @seealso \code{\link[BAMBI]{dvmsin}}

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
