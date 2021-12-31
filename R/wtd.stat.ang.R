#' Weighted extrinsic mean direction and mean resultant length
#'
#' \code{wtd.stat.ang} computes weighted extrinsic mean direction and
#'   mean resultant length.
#'
#' @param data angular data whose elements are in \eqn{[0, 2\pi)}
#' @param w numeric vector whose each element is non-negative and
#'   \code{sum(w) == 1}. Moreover, the length of \code{w} is the same with
#'   \code{nrow(data)}.
#' @return list which is consisting of the following components:
#'
#'   \code{Mean} weighted extrinsic mean direction
#'
#'   \code{R} mean resultant length
#' @export
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
#'
#' @examples
#' data <- matrix(c(pi/3, pi/3, pi/2,
#'                  pi, pi/4, pi/2,
#'                  0, pi/3, pi/6),
#'                ncol = 3, byrow = TRUE)
#' w <- c(0.3, 0.3, 0.4)
#' wtd.stat.ang(data, w)
wtd.stat.ang <- function(data, w){
  # computes weighted extrinsic mean direction and mean resultant length

  if(is.vector(data)){ data <- t(as.matrix(data))}
  # note w is multiplied to each column of the former
  wtd_ext_mean <- colMeans((cos(data) + 1i * sin(data)) * w)

  ang <- Arg(wtd_ext_mean)
  ang <- ifelse(ang < 0, ang + 2 * pi, ang)
  R <- Mod(wtd_ext_mean)

  return(list(Mean = ang, R = R))
}
