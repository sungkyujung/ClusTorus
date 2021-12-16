#' K-Means Clustering to K-Spheres Clustering on Torus
#'
#' \code{ellip.kmeans.torus} prepares the parameters for conformity scores
#'   which are derived by k-means clustering on torus.
#'
#' @param data data n x d matrix of toroidal data on \eqn{[0, 2\pi)^d}
#' @param centers either the number of clusters or a set of initial
#'   cluster centers. If a number, a random set of row in x is
#'   chosen as the initial centers.
#' @param type character which must be "homogeneous-circular",
#'  "heterogeneous-circular", or "general".
#'   If "homogeneous-circular", the radii of k-spheres are identical.
#'   If "heterogeneous-circular", the radii of k-spheres may be different.
#'   If "ellipsoids", cluster with k-ellipsoids without optimized parameters.
#'   If, "general", clustering with k-ellipsoids. The parameters to construct
#'   the ellipses are optimized with elliptical k-means algorithm, which is
#'   modified for toroidal space. See references for the detail.
#'   Default is "homogeneous-circular".
#' @param init determine the initial parameter for option "general". Must be
#'   "kmeans" or "hierarchical".
#'   If "kmeans", the initial parameters are obtained with extrinsic kmeans
#'   method.
#'   If "hierarchical", the initial parameters are obtained with hierarchical
#'   clustering method. Default is "hierarchical".
#' @param ... Further arguments for argument \code{init}. If \code{init = "kmeans"},
#'   these are for \code{\link[stats]{kmeans}}. If \code{init = "hierarchical"},
#'   there are for \code{\link[stats]{hclust}}.
#' @param additional.condition boolean index.
#'   If \code{TRUE}, a singular matrix will be altered to the scalar identity.
#' @param THRESHOLD number of threshold for difference between updating and
#'   updated parameters. Default is 1e-10.
#' @param maxiter the maximal number of iteration. Default is 200.
#' @param verbose boolean index, which indicates whether display
#'   additional details as to what the algorithm is doing or
#'   how many loops are done. Default is \code{TRUE}.
#'
#' @return returns a list,
#'   containing all values which determines the shape and
#'   location of spheres.
#' @export
#' @seealso
#'   \code{\link[ClusTorus]{kmeans.torus}}
#' @references S. Jung, K. Park, and B. Kim (2021),
#'   "Clustering on the torus by conformal prediction", and
#'   Jaehyeok Shin, Alessandro Rinaldo and Larry Wasserman (2019),
#'   "Predictive Clustering"
#' @examples
#' data <- ILE[1:200, 1:2]
#'
#' ellip.kmeans.torus(data, centers = 3, type = "general", init = "hierarchical")
ellip.kmeans.torus <- function(data, centers = 10,
                            type = c("homogeneous-circular",
                                     "heterogeneous-circular",
                                     "ellipsoids",
                                     "general"),
                            init = c("kmeans", "hierarchical"),
                            additional.condition = TRUE,
                            THRESHOLD = 1e-10, maxiter = 200,
                            verbose = TRUE, ...){

  # Returns a sphere.param object, containing all values which determines
  # the shape and location of spheres

  # type determines kmeans-fitting method. If "identical", the radii of
  # shperes are the same, and if not, the radii may be different.
  if (is.null(type)){ type <- "homogeneous-circular" }
  if (is.null(init)){ init <- "hierarchical" }

  type <- match.arg(type)
  init <- match.arg(init)
  d <- ncol(data)
  n <- nrow(data)

  sphere.param <- list(mu = NULL, Sigmainv = NULL, c = NULL)

  # Use extrinsic kmeans clustering for initial center points.
  # centers is given as a number, in default, but it may also be given
  # as a matrix which indicates the toroidal points.

  # -------------- initializing ----------------
  if (init == "kmeans"){
    kmeans.out <- kmeans.torus(data, centers, ...)
  } else {
    J <- ifelse(is.null(ncol(centers)), centers, ncol(centers))
    kmeans.out <- hcluster.torus(data, J = centers, ...)
  }

  centroid <- kmeans.out$centers
  J <- nrow(centroid)

  # 1. homogeneous spheres
  # In fact, the initialized parameters are for the identical case.
  # sphere.param$mu1 <- centroid[, 1]
  # sphere.param$mu2 <- centroid[, 2]
  sphere.param$mu <- centroid
  sphere.param$c <- rep(0, J)

  for(j in 1:J){
    sphere.param$Sigmainv[[j]] <- diag(d)
  }


  # 2. heterogeneous spheres --------------------------
  if (type == "heterogeneous-circular"){

    for(j in 1:J){

      # if the size of cluster is 1, the cluster contains only one point.
      nj <- kmeans.out$size[j]
      pi_j <- nj / n
      sigma_j <- ifelse(kmeans.out$size[j] <= 1,
                        1e-6, kmeans.out$withinss[j] / (nj * d))

      sphere.param$c[j] <- 2 * log(pi_j) - d * log(sigma_j)
      sphere.param$Sigmainv[[j]] <- diag(d) / sigma_j
    }
  }

  # 3. kmeans to ellipsoids ----------------------------
  else if (type == "ellipsoids") {

    for (j in 1:J){

      nj <- kmeans.out$size[j]
      pi_j <- nj / n

      dat.j <- data[kmeans.out$membership == j, ]

      # z <- tor.minus(dat.j, c(sphere.param$mu1[j], sphere.param$mu2[j]))
      z <- tor.minus(dat.j, sphere.param$mu[j, ])

      S <- t(z) %*% z / nrow(z)

      # additional assumption to S : axis-aligned
      if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
        S <- diag(diag(S))
      }

      # additional assumption to S : sphere
      # only implemented when verbose == TRUE
      if (additional.condition){
        if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
          cnt.singular <- cnt.singular + 1
          S <- sum(S) / d * diag(d)
        }
      }

      # vanishing the ellipsoid even if the additional condition is given.
      if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
        S <- 1e-6 * diag(d)
      }

      sphere.param$Sigmainv[[j]] <- solve(S)

      # Step.5 -----------------------------------
      pi_j <- ifelse(sum(kmeans.out$membership == j) == 0,
                     1e-6, sum(kmeans.out$membership == j) / n)
      # update c's
      sphere.param$c[j] <- 2 * log(pi_j) - log(det(S))
    }
  }


  # 4. general ellipsoids ------------------------
  # initialize the parameters with EMsinvMmix.init and norm.appr.param
  # Use generalized Lloyd's algorithm
  else if (type == "general"){
    # Step.1 --------------------------------------------
    # initialize the parameters


    for (j in 1:J){

      nj <- kmeans.out$size[j]
      pi_j <- nj / n

      # dat.j <- data[kmeans.out$membership == j, ]
      # z <- tor.minus(dat.j, c(sphere.param$mu1[j], sphere.param$mu2[j]))
      z <- tor.minus(data[kmeans.out$membership == j, ], sphere.param$mu[j, ])

      S <- t(z) %*% z / nrow(z)

      # additional assumption to S : axis-aligned
      if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
        S <- diag(diag(S))
      }

      # additional assumption to S : sphere
      # only implemented when verbose == TRUE
      if (additional.condition){
        if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
          S <- sum(S) / d * diag(d)
        }
      }

      # vanishing the ellipsoid even if the additional condition is given.
      if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
        S <- 1e-6 * diag(d)
      }

      sphere.param$Sigmainv[[j]] <- solve(S)

      # Step.5 -----------------------------------
      pi_j <- ifelse(sum(kmeans.out$membership == j) == 0,
                     1e-6, sum(kmeans.out$membership == j) / n)

      # update c's
      sphere.param$c[j] <- 2 * log(pi_j) - log(det(S))

    }

    # vectorize the sphere.param: this will be used for escaping loop
    param.seq <- unlist(sphere.param)

    if (verbose){
      cat("ellip.kmeans.torus: fitting parameters with option ",type, ", J =", J, "\n")
    }

    cnt <- 1
    wmat <- ehatj <- matrix(0, n, J)

    while(TRUE){
      cnt <- cnt + 1

      if(verbose){if (cnt %% 5 == 0){cat(".")}}

      # Step.2 ------------------------------------
      # prepare w's which work like weights

      ehatj <- ehat.eval(data, sphere.param)
      # maxj <- apply(ehatj, 1, which.max)
      # evaluate wmat
      for(j in 1:J){ wmat[, j] <- max.col(ehatj, ties.method = "first") == j }

      # Step.3 -------------------------------------
      # update mu's

      # wmat.mul <- apply(wmat, 2, '*', data)
      # wmat.mul <- rbind(wmat.mul, wmat)
      # wmat.mul <- apply(wmat.mul, 2, function(x){
      #   wtd.stat.ang(matrix(x[1:(d * n)], n, byrow = F),
      #                w = x[((d * n) + 1):length(x)] / sum(x[((d * n) + 1):length(x)]))$Mean})
      wmat.mul <- apply(wmat, 2, function(x){
        dat.j <- data[x == 1, ]
        nj <- length(dat.j) / d
        if(nj > 0){
          return(wtd.stat.ang(dat.j, w = rep(1, nj) / nj)$Mean)
        } else { return(rep(0, d)) }
      })


      # sphere.param$mu1 <- mu[, 1]
      # sphere.param$mu2 <- mu[, 2]
      sphere.param$mu <- t(wmat.mul)

      # Step.4 and Step.5

      for (j in 1:J){
        # dat.j <- data[wmat[, j] == 1, ]
        # Step.4 -----------------------------------
        # z <- tor.minus(dat.j, c(sphere.param$mu1[j], sphere.param$mu2[j]))
        z <- tor.minus(data[wmat[, j] == 1, ], sphere.param$mu[j, ])

        # evaluate the MLE of Sigma_j
        S <- t(z) %*% z / nrow(z)

        # additional assumption to S : axis-aligned
        if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
          S <- diag(diag(S))
        }

        # additional assumption to S : sphere
        # only implemented when additional.condition == TRUE
        if (additional.condition){
          if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
            S <- sum(S) / d * diag(d)
          }
        }

        # vanishing the ellipsoid even if the additional condition is given.
        if (det(S) < THRESHOLD || sum(is.na(S)) != 0){
          S <- 1e-6 * diag(d)
        }

        sphere.param$Sigmainv[[j]] <- solve(S)

        # Step.5 -----------------------------------
        pi_j <- ifelse(sum(wmat[, j]) == 0, 1e-6, sum(wmat[, j]) / n)

        # update c's
        sphere.param$c[j] <- 2 * log(pi_j) - log(det(S))
      }

      diff <- sum((param.seq - unlist(sphere.param))^2, na.rm = TRUE)
      param.seq <- unlist(sphere.param)

      if (cnt >= maxiter | diff < THRESHOLD){
        if(verbose){
          cat("Done")
          cat("\n")
        }
        break}
    }
    sphere.param$loglkhd <- 0.5 * sum(do.call(pmax,
                                               as.data.frame(ehat.eval(data, sphere.param)))) - n * d * log(2 * pi) / 2

    sphere.param$singular <- c()
    for (j in 1:J){
      if (sphere.param$Sigmainv[[j]][1, 1] == 1e+6){
        sphere.param$singular <- c(sphere.param$singular, j)
      }
    }
  }
  return(sphere.param)
}
