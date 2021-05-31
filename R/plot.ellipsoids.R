# Plot Ellipsoids on 2-dimensional Toroidal Space
#
# \code{ploting.ellipsoids} generates a 2-dimensional plot for the given ellipsoids,
#   with prespecified coordinates.
#
# @param data data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
# @param ellipse.param list which is consisting of mean of each angular
#   coordinate, inverse of each covariance matrix, and constant term
# @param t a numeric value which determines the size of ellipses.
# @param coord a 2-vector for prespecifing the coordinates.
#   Default value is c(1, 2).
# @return a ggplot2-based 2-dimensional plot with ellipses
# @seealso \code{\link[ClusTorus]{icp.torus.score}}, \code{\link[ClusTorus]{cluster.assign.torus}}
# @examples
# \dontrun{
# parammat <- matrix(c(0.4, 0.3, 0.3,
#                      20, 25, 25,
#                      30, 25, 20,
#                      1, 2, 3,
#                      1, 2, 3,
#                      0, 2, 4), nrow = 6, byrow =TRUE)
#
# ellipse.param <- norm.appr.param(parammat)
#
# coord <- c(1, 3)
# t <- 0.5
#
# plot.ellipsoids(ellipse.param, t, coord)
# }

ploting.ellipsoids <- function(data, ellipse.param, t, coord = c(1, 2)){
  if (is.vector(coord)) {coord <- t(as.matrix(coord))}
  if (ncol(coord) != 2 | !is.numeric(coord))
  {stop("Invalid coordinates: coord must be a n x 2-dimensional numeric vector/matrix.")}

  d <- ncol(data)
  angle.names <- colnames(data)
  if (max(coord) > d | min(coord) < 1)
  {stop("Invalid coordinates: Out of dimensionality.")}

  plot_list <- list()

  J <- length(ellipse.param$c)
  theta <- seq(0, 2 * pi, length.out = 199)
  Z <- cbind(cos(theta), sin(theta))

  shift <- matrix(0, ncol = 2, nrow = 9)
  shift[, 1] <- c(0, 2 * pi, -2 * pi)
  shift[, 2] <- rep(c(0, 2 * pi, -2 * pi), each = 3)

  for (i in 1:nrow(coord)){
    plot_list[[i]] <- local({
      i <- i
      ang1 <- coord[i, 1]
      ang2 <- coord[i, 2]
      g2 <- ggplot2::ggplot()

      for (j in 1:J){
        mu <- ellipse.param$mu[j, c(ang1, ang2)]
        Sinv <- ellipse.param$Sigmainv[[j]][c(ang1, ang2), c(ang1, ang2)]
        c.minus.t <- ellipse.param$c[j] - t

        if (c.minus.t <= 0) {next}

        M <- eigen(Sinv/c.minus.t)
        Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
        R <- Mmhalf %*% t(Z)
        for( shift.id in 1:9){
          RR <- R + mu + shift[shift.id,]
          plot.data <- data.frame(angle1 = RR[1,], angle2 = RR[2,], value = 1)
          g2 <- g2 + ggplot2::geom_polygon(ggplot2::aes(x = .data$angle1, y = .data$angle2),
                                           color = "blue",alpha = 0.1,
                                           data = plot.data) +
            ggplot2::xlab(angle.names[ang1]) + ggplot2::ylab(angle.names[ang2])
        }
      }
      g2 <- g2 + ggplot2::scale_x_continuous(breaks = c(0,1,2)*pi,
                           labels = c("0","pi", "2pi")) +
        ggplot2::scale_y_continuous(breaks = c(0,1,2)*pi,
                           labels = c("0","pi","2pi")) +
        ggplot2::coord_cartesian(xlim = c(0, 2*pi), ylim = c(0, 2*pi), expand = FALSE) +
        ggplot2::geom_point(ggplot2::aes(x = data[, ang1], y = data[, ang2]))
    })
  }

  plot_list
}
