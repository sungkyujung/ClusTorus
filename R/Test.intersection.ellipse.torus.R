# Intersection of two ellipses on torus
#
# \code{Test.intersection.ellipse.torus} evaluates whether two ellipses
#   on torus intersect.
#
# @param ellipse.param list which is consisting of mean of each angular
#   coordinate, inverse of each covariance matrix, and constant terms
# @param index 2-dimensional vector which indicates the ellipses that
#   we will check.
# @param t a numeric value which determines the size of ellipses.
# @return If they intersect, then return \code{TRUE}. If not,
#   then return \code{FALSE}.
# @seealso \code{\link[ClusTorus]{Test.intersection.ellipse}}
# @references S. Jung, K. Park, B. Kim (2020), "Clustering on the torus
#   by conformal prediction"
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
# index <- c(1, 3)
# t <- 0.5
#
# Test.intersection.ellipse.torus(ellipse.param, index, t)
# }
Test.intersection.ellipse.torus <- function(ellipse.param, index, t){

  i <- index[1]
  j <- index[2]

  d <- ncol(ellipse.param$Sigmainv[[1]])

  # mean.1 <- matrix(c(ellipse.param$mu1[i], ellipse.param$mu2[i]),ncol = 1)
  mean.1 <- ellipse.param$mu[i, ]
  Sinv1 <- ellipse.param$Sigmainv[[i]]
  c1.minus.t <- ellipse.param$c[i] - t

  # mean.2 <- matrix(c(ellipse.param$mu1[j], ellipse.param$mu2[j]),ncol = 1)
  mean.2 <- ellipse.param$mu[j, ]
  Sinv2 <- ellipse.param$Sigmainv[[j]]
  c2.minus.t <- ellipse.param$c[j] - t

  if(c1.minus.t <= 0 || c2.minus.t <= 0){
    overlap <- FALSE
    return(overlap)
  }

  M.1 <- Sinv1 / c1.minus.t
  M.2 <- Sinv2 / c2.minus.t

  # preparing for case 1
  a <- M.1[1, 1]
  b <- M.2[1, 1]

  # preparing for case 2

  r1 <- 1 / eigen(M.1)$values[d]
  r2 <- 1 / eigen(M.2)$values[d]

  # generate shifting vectors
  shift <- matrix(0, ncol = d, nrow = 3^d)
  for (i in 1:d){
    shift[, i] <- rep(c(0, 2 * pi, -2 * pi), each = 3^(i - 1))
  }


  # case 1 : both ellipsoids are spheres
  if ((sum(M.1/a) == d) && (sum(M.2/b) == d)){

    center.dist <- sqrt(sum(ang.minus(mean.1, mean.2)^2))
    radius.sum <- sqrt(1/a) + sqrt(1/b)

    overlap <- center.dist <= radius.sum

    # case 2 : sum of the longest radii of ellipsoids are smaller than pi
  } else if (sqrt(r1) + sqrt(r2) <= pi){

    mu <- as.vector(ang.minus(mean.2, mean.1))
    dist <- rowSums(sweep(shift, 2, mu)^2)
    shift.ind <- which.min(dist)

    overlap <- Test.intersection.ellipse(shift[shift.ind, ], M.1, mu, M.2)

    # mean.1.shift <- sweep(shift, 2, mean.1, "+")
    # shift.dist <- rowSums(sweep(mean.1.shift, 2, mean.2)^2)
    # shift.ind <- which.min(shift.dist)
    #
    # overlap <- Test.intersection.ellipse(mean.1.shift[shift.ind], M.1, mean.2, M.2)

    # general case
  } else {

    # overlap <- FALSE
    # # method 1 : using for loop ----------------
    # for(trial in 1:3^d){
    #   overlap <- Test.intersection.ellipse(mean.1, M.1, mean.2 + shift[trial,], M.2)
    #   if (overlap) {break}
    # }

    # return(overlap)

    # method 2 : using Vectorize ----------------
    # Test <- Vectorize(function(i) {Test.intersection.ellipse(mean.1, M.1, mean.2 + shift[i, ], M.2)})
    #
    # overlap <- sum(Test(1:3^d))
    # return(overlap >= 1)

    # method 3 : using purrr::map ---------------

    overlap.results <- purrr::map_int(1:dim(shift)[1], function(i)
      {Test.intersection.ellipse(mean.1, M.1, mean.2 + shift[i, ], M.2)})

    overlap <- sum(overlap.results) >= 1
  }

  return(overlap)
}
