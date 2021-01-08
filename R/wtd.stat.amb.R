# Weighted means in the ambient space and the cross matrix
#
# \code{wtd.stat.amb} returns the vector-valued weighted means in the
#   ambient space and the cross matrix.
#
# @description See section A.1 of the article 'Clustering on the
#   torus by conformal prediction', S. Jung, K. Park, and B. Kim (2020)
#
# @inheritParams wtd.stat.ang
# @return list which is consisting of \code{y1bar}, \code{y2bar}, \code{S12}.
# @references 'S. Jung, K. Park, and B. Kim (2020),
#   "Clustering on the torus by conformal prediction"
#

wtd.stat.amb <- function(data, w){
  # returns the vector-valued weighted means in the ambient space and the cross matrix

  # note w is multiplied to each column of the former

  y <- cos(data) + 1i * sin(data)
  wtd_ext_mean <- colMeans( y * w )


  ymean <- rbind( Re(wtd_ext_mean),
                  Im(wtd_ext_mean))
  #rownames(ymean) <- c("coord.1","coord.2")
  #colnames(ymean) <- c("mu1","mu2")

  y1bar <- ymean[,1]
  y2bar <- ymean[,2]

  S12 <-    colMeans(
    cbind(Im(y[,1]) * Im(y[,2]) * w,
          -Re(y[,1]) * Im(y[,2]) * w,
          -Im(y[,1]) * Re(y[,2]) * w,
          Re(y[,1]) * Re(y[,2]) * w) )
  S12 <- matrix(S12, nrow = 2)

  return( list(y1bar = y1bar, y2bar = y2bar, S12 = S12))
}
