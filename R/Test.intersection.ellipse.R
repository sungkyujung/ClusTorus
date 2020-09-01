#' Intersection of two ellipses on R^2
#'
#' \code{Test.intersection.ellipse} evaluates whether two ellipses
#'   on R^2 intersect.
#'
#' @param mean.1,mean.2 2-dimensional vectors which indicate the
#'   centers of each ellipse.
#' @param M.1,M.2 2 x 2 matrices which determine the shape of
#'   each ellipse.
#' @return If they intersect, then return \code{TRUE}. If not,
#'   then return \code{FALSE}.
#' @export
#' @references We use the efficient algorithm given by
#'  David Eberly (2019), "Intersection of Ellipses",
#'  Geometric Tools, Redmond WA 98052
#' @seealso \code{\link[ClusTorus]{Test.intersection.ellipse.torus}}
#' @export
#' @examples
#' \dontrun{
#' mean.1 <- c(1, 0)
#' mean.2 <- c(0, 2)
#' M.1 <- diag(c(3, 4))
#' M/2 <- matrix(c(5, 2, 2, 1), nrow = 2, byrow = TRUE)
#'
#' Test.intersection.ellipse(mean.1, M.1, mean.2, M.2)
#' }
Test.intersection.ellipse <- function(mean.1, M.1, mean.2, M.2){
  #
  #
  #   i <- index[1]
  #   j <- index[2]
  #
  #   mean.1 <- matrix(c(ellipse.param$mu1[i], ellipse.param$mu2[i]),ncol = 1)
  #   Sinv1 <- ellipse.param$Sigmainv[[i]]
  #   c1.minus.t <- ellipse.param$c[i] - conf.score.level
  #
  #   mean.2 <- matrix(c(ellipse.param$mu1[j], ellipse.param$mu2[j]),ncol = 1)
  #   Sinv2 <- ellipse.param$Sigmainv[[j]]
  #   c2.minus.t <- ellipse.param$c[j] - conf.score.level
  #
  #   if(c1.minus.t <= 0 || c2.minus.t <= 0){
  #     Ind.Overlap <- 0
  #     return(Ind.Overlap)
  #   }
  #
  #   M.1 <- Sinv1 / c1.minus.t
  #   M.2 <- Sinv2 / c2.minus.t

  # Now, each ellipse satisfies (x - mean.j)^T M.j (x - mean.j) = 1

  Eig.1 <- eigen(M.1)
  R.1 <- Eig.1$vectors
  D.1 <- diag((Eig.1$values))
  Eig.2 <- eigen(M.2)
  R.2 <- Eig.2$vectors
  D.2 <- diag((Eig.2$values))

  # Reduction to the unit circle and axis-aligned ellipse
  K.3 <- sqrt(D.1) %*% t(R.1) %*% (mean.2 - mean.1)
  M.3 <- diag((1/sqrt(diag(D.1)))) %*% t(R.1) %*% R.2 %*% D.2 %*% t(R.2) %*% R.1 %*% diag((1/sqrt(diag(D.1))))

  Eig.3 <- eigen(M.3)
  R <- Eig.3$vectors[,order(Eig.3$values)]
  D <- diag(Eig.3$values[order(Eig.3$values)])
  K <- t(R) %*% K.3

  # package 'polynom' is required
  d0 <- D[1,1]
  d1 <- D[2,2]
  k0 <- K[1]
  k1 <- K[2]
  coef.0 <- 1 - d0 * k0^2 - d1 * k1^2
  coef.1 <- -2 * d0 - 2 * d1 + 2 * d0 * d1 * k0^2 + 2 * d0 * d1 * k1^2
  coef.2 <- d0^2 + d1^2 + 4*d0*d1 - d0*d1^2*k0^2 - d0^2*d1*k1^2
  coef.3 <- -2 * d0 * d1^2 - 2 * d0^2 * d1
  coef.4 <- d0^2 * d1^2
  f.s <- polynom::polynomial(c(coef.0, coef.1, coef.2, coef.3, coef.4))
  root.f.s <- Re(solve(f.s))

  P.0 <- c(d0 * k0 * min(root.f.s) / (d0 * min(root.f.s) - 1), d1 * k1 * min(root.f.s) / (d1 * min(root.f.s) - 1))
  P.1 <- c(d0 * k0 * max(root.f.s) / (d0 * max(root.f.s) - 1), d1 * k1 * max(root.f.s) / (d1 * max(root.f.s) - 1))

  # From now on, test whether two ellipses are overlapped or not
  minDistance <- sqrt(sum(P.0^2))
  maxDistance <- sqrt(sum(P.1^2))
  Ind.Overlap <- 0
  if (maxDistance <= 1){
    Ind.Overlap <- 1
  }else{
    if (minDistance < 1){
      Ind.Overlap <- 1
    }else if (minDistance > 1){
      if (d0 * k0^2 + d1 * k1^2 - 1 > 0){
        Ind.Overlap <- 0
      }else{
        Ind.Overlap <- 1
      }
    }else{
      if (d0 * k0^2 + d1 * k1^2 - 1 > 0){
        Ind.Overlap <- 0
      }else{
        Ind.Overlap <- 1
      }
    }
  }
  return(Ind.Overlap == 1)
}
