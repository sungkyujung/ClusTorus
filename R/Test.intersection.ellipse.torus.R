#' Intersection of two ellipses on torus
#'
#' \code{Test.intersection.ellipse.torus} evaluates whether two ellipses
#'   on torus intersect.
#'
#' @param ellipse.param list which is consisting of mean of each angular
#'   coordinate, inverse of each covariance matrix, and constant terms
#' @param index 2-dimensional vector which indicates the ellipses that
#'   we will check.
#' @param t a numeric value which determines the size of ellipses.
#' @return If they intersect, then return \code{TRUE}. If not,
#'   then return \code{FALSE}.
#' @seealso \code{\link[ClusTorus]{Test.intersection.ellipse}}
#' @references S. Jung, K. Park, B. Kim (2020), "Clustering on the torus
#'   by conformal prediction"
#' @export
#' @examples
#' \dontrun{
#' parammat <- matrix(c(0.4, 0.3, 0.3,
#'                      20, 25, 25,
#'                      30, 25, 20,
#'                      1, 2, 3,
#'                      1, 2, 3,
#'                      0, 2, 4), nrow = 6, byrow =TRUE)
#'
#' elipse.param <- norm.appr.param(parammat)
#'
#' index <- c(1, 3)
#' t <- 0.5
#'
#' Test.intersection.ellipse.torus(ellipse.param, index, t)
#' }
Test.intersection.ellipse.torus <- function(ellipse.param, index, t){

  i <- index[1]
  j <- index[2]

  mean.1 <- matrix(c(ellipse.param$mu1[i], ellipse.param$mu2[i]),ncol = 1)
  Sinv1 <- ellipse.param$Sigmainv[[i]]
  c1.minus.t <- ellipse.param$c[i] - t

  mean.2 <- matrix(c(ellipse.param$mu1[j], ellipse.param$mu2[j]),ncol = 1)
  Sinv2 <- ellipse.param$Sigmainv[[j]]
  c2.minus.t <- ellipse.param$c[j] - t

  if(c1.minus.t <= 0 || c2.minus.t <= 0){
    overlap <- FALSE
    return(overlap)
  }

  M.1 <- Sinv1 / c1.minus.t
  M.2 <- Sinv2 / c2.minus.t


  shift <- matrix(0,ncol = 2, nrow = 9)
  shift[,1] <- c(0,2*pi,-2*pi)
  shift[,2] <- rep(c(0,2*pi,-2*pi), each = 3)


  shift.id <- 1
  overlap <- FALSE
  for(trials in 1:9){
    overlap <- Test.intersection.ellipse(mean.1, M.1, mean.2 + shift[shift.id,], M.2)
    ifelse(overlap, break, shift.id <- shift.id +1)
  }
  return(overlap)

}