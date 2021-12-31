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
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
#'
#'   Di Marzio, M., Panzera, A., & Taylor, C. C. (2011). Kernel density estimation on the torus. \emph{Journal of Statistical Planning and Inference}, 141(6), 2156-2173.
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
  # if data contains NAs, the rows containing NAs are removed by na.omit()

  data <- stats::na.omit(data)

  if (!is.matrix(data)) {data <- as.matrix(data)}

  N <- nrow(eval.point)
  n <- nrow(data)

  eval.point.bind <- rbind(eval.point, data)

  concentration = concentration[1]
  if (!is.numeric(concentration) | concentration <= 0) {
    concentration <- 25
    warning("Concentration must be a positive number. Reset as concentration = 25 (default)\n")
  }

  phat <- kde.torus(data, eval.point.bind,
                    concentration = concentration)

  phat.grid <- phat[1:N]                 # hat{p}(eval.point)
  phat.data <- sort( phat[(N + 1):(N + n)], index.return = TRUE)
  data <- data[phat.data$ix, ] # data reordered to satisfy y_i = y_(i)
  phat.data <- phat.data$x    # hat{p}(y_(i)) sorted (increasing)

  cp.torus <- NULL


  if (!is.null(level)){

    if (level < 0 || level > 1) {
      level <- 0.1
      warning("Level must be numeric and between 0 and 1. Reset as level = 0.1 (default)")
    }

    nalpha <- length(level)

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

  structure(list(concentration = concentration, level = level,
                 cp.torus = cp.torus,
                 grid = eval.point,
                 phat.grid = phat.grid,
                 phat.data = phat.data,
                 data.sorted = data), class = "cp.torus.kde")
}
