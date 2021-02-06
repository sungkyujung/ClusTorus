#' Conformal prediction set indices with kernel density estimation
#'
#' \code{cp.torus.kde} computes conformal prediction set indices
#'   (TRUE if in the set) using kernel density estimation as conformity score.
#'
#' @inheritParams kde.torus
#' @param level either a scalar or a vector, or even \code{NULL}. Default value
#'   is 0.1.
#' @return If \code{level} is \code{NULL}, then return kde at \code{eval.point}
#'   and at data points.
#'
#'   If \code{level} is a vector, return the above and prediction set indices
#'   for each value of level.
#' @seealso \code{\link{kde.torus}}, \code{\link{grid.torus}}
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction"
#' @export
#' @examples
#' data <- ILE[1:200, 1:2]
#' cp.torus.kde(data, eval.point = grid.torus(),
#'              level = 0.05, concentration = 25)

cp.torus.kde <- function(data, eval.point = grid.torus(),
                         level= 0.1,
                         concentration = 25){
  # Computes conformal prediction set indices (TRUE if in the set) using kde as conformal scores
  # level can be either a scalar or a vector, or even null.
  # if level is null, return kde at eval.point and at data points.
  # if level is a vector, return the above and Prediction set indices for each value of level.
  if (!is.matrix(data)) {data <- as.matrix(data)}

  N <- nrow(eval.point)
  n <- nrow(data)

  eval.point.bind <- rbind(eval.point, data)
  #
  phat <- kde.torus(data, eval.point.bind,
                    concentration = concentration)

  phat.grid <- phat[1:N]                 # hat{p}(eval.point)
  phat.data <- sort( phat[(N + 1):(N + n)], index.return = TRUE)
  data <- data[phat.data$ix, ] # data reordered to satisfy y_i = y_(i)
  phat.data <- phat.data$x    # hat{p}(y_(i)) sorted (increasing)

  nalpha <- length(level)
  cp.torus <- NULL

  if (!is.null(level)){
    for (i in 1:nalpha){
      ialpha <- floor( (n + 1) * level[i])
      # indices for inclusion in L-
      Lminus <- phat.grid >= phat.data[ialpha]

      # indices for inclusion in L+
      const_vM2 <- (2 * pi * besselI(concentration, 0))^2
      zeta <- (exp(concentration * 2) - exp(-concentration * 2)) / const_vM2
      Lplus <- phat.grid >=  phat.data[ialpha] - zeta / n

      # indices for inclusion in Cn(alpha)
      K <- (exp(concentration * 2) -
              exp(concentration * (  rowSums(cos(sweep(eval.point, 2, data[ialpha, ], "-"))  )))) /
        const_vM2 # K(0) - K(y_(ialpha) - y)
      Cn <- phat.grid - phat.data[ialpha] + K/n >= 0

      cp.torus.i <- data.frame(phi = eval.point[, 1], psi = eval.point[, 2],
                               Lminus = Lminus, Cn = Cn, Lplus = Lplus, level = level[i])
      cp.torus <- rbind(cp.torus, cp.torus.i)
    }
  }


  list(cp.torus = cp.torus,
       grid = eval.point,
       phat.grid = phat.grid,
       phat.data = phat.data,
       data.sorted = data
  )
}
