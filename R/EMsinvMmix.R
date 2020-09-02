EMsinvMmix <- function(data, J = 4, parammat = EMsinvMmix.init(data, J),
                       THRESHOLD = 1e-10, maxiter = 200,
                       type = c("circular","axis-aligned","general"),
                       kmax = 500,
                       verbose = TRUE){
  # input data: n x 2 angles
  # input parammat: 6 x J matrix of initial values
  # input J: number of components
  # R is the number of seq. updates within M-step.
  # return fitted parameters of J-mixture of bivariate sine von Mises

  # initialize
  # maxiter <- 100
  # THRESHOLD <- 1e-10
  type <- match.arg(type)
  n <- nrow(data)
  param.seq <- as.vector(parammat)

  loglkhd.seq <- sum( BAMBI::dvmsinmix(data,kappa1 = parammat[2,],
                                       kappa2 = parammat[3,],
                                       kappa3 = parammat[4,],
                                       mu1 = parammat[5,],
                                       mu2 = parammat[6,],
                                       pmix = parammat[1,],log = TRUE))

  pimat <- matrix(0,nrow = n, ncol = J)
  #cat(ncol(parammat))
  if(verbose){cat("EMsinvMmix: fitting vM2 with option ",type, ", J=", J, ".\n")}

  cnt <- 1
  while(TRUE){
    cnt <- cnt + 1

    if(verbose){if (cnt %% 10 == 0){ cat(cnt)}}

    # E-step ------------------------------------------------------------------
    # E-step to update pimat = Prob of ith sample to be in the jthe group
    for (j in 1:J){
      # for each component,
      # compute the n-vector of conditional density: jth column of pimat
      pimat[,j] <- BAMBI::dvmsin(data,kappa1 = parammat[2,j],
                                 kappa2 = parammat[3,j],
                                 kappa3 = parammat[4,j],
                                 mu1 = parammat[5,j],
                                 mu2 = parammat[6,j]
      ) * parammat[1,j]
    }
    # expected class membership:
    pimat <- ( pimat  / rowSums(pimat) ) # for each sample, i (each row), sum to 1
    # rowSums(pimat): all ones.


    # M-step ------------------------------------------------------------------
    # Given pimat, update parammat.
    # For pi (group probability), update at once first
    # For kappas and mus, separately update for each component, depending on the type of model.
    #

    # M-step 1. pi_j (j=1,...,J) update
    parammat[1,] <- colMeans(pimat)

    # M-step 2. vM parameters update
    if ( type == "axis-aligned"){
      # This is case I where lambda = 0

      # For each j, compute weighted angluar mean and mean resultant length
      for(j in 1:J){
        pivec <- pimat[,j]
        wstat <-wtd.stat.ang(data, w = pivec)
        parammat[2,j] <- min( c( Ainv(wstat$R[1] / parammat[1,j]), kmax))
        parammat[3,j] <- min( c( Ainv(wstat$R[2] / parammat[1,j]), kmax))
        parammat[4,j] <- 0
        parammat[5:6,j] <- wstat$Mean
      }

    } else if (type == "circular"){
      # This is case II where lambda = 0 and kappa1= kappa2= kappa.

      # For each j, compute weighted angluar mean and a joint mean resultant length
      for(j in 1:J){
        pivec <- pimat[,j]
        wstat <-wtd.stat.ang(data, w = pivec)
        kappa_next <- min( c( Ainv( mean(wstat$R) / parammat[1,j] )  , kmax))
        parammat[2,j] <- kappa_next
        parammat[3,j] <- kappa_next
        parammat[4,j] <- 0
        parammat[5:6,j] <- wstat$Mean
      }



    } else {
      # This is the general case. Apply Conditional Maximization Either.

      # Step (a) is done (Step 1)

      # Step (b) Update mean parameter
      for (j in 1:J){
        # initialize
        mu1 <- c(cos(parammat[5,j]),sin(parammat[5,j]))
        mu2 <- c(cos(parammat[6,j]),sin(parammat[6,j]))

        # compute statistics
        pivec <- pimat[,j]
        wstat <- wtd.stat.amb(data, w = pivec)

        # repeatedly update using ybar, S12, kappa1, kappa2, lambda, and mu1,mu2
        muupdate <- sinvM.ECMEb(wstat,
                                kappa1 = parammat[2,j],
                                kappa2 = parammat[3,j],
                                lambda = parammat[4,j],
                                mu1 = parammat[5,j],
                                mu2 = parammat[6,j])
        parammat[5,j] <- muupdate$mu1
        parammat[6,j] <- muupdate$mu2
      }

      # Step (c) Update concentration kappa1 kappa2

      # by maximizing the observed data log-likelihood

      initialkappas <- c(parammat[2,],parammat[3,])
      a<- optim(par =  initialkappas,
                fn = function(kappas){
                  sum( BAMBI::dvmsinmix(data,kappa1 = kappas[1:J],
                                        kappa2 = kappas[(J+1):(2*J)],
                                        kappa3 = parammat[4,],
                                        mu1 = parammat[5,],
                                        mu2 = parammat[6,],
                                        pmix = parammat[1,],log = TRUE))
                },
                lower = rep(1e-10,2*J),
                upper = rep(kmax,2*J),
                method = "L-BFGS-B",
                control = list(fnscale = -1))
      parammat[2,] <- a$par[1:J]
      parammat[3,] <- a$par[(J+1):(2*J)]

      # Step (d) Update lambdas
      # compute the bounds for lambdas
      lbound <- sqrt(parammat[2,]*parammat[3,])
      a <- optim(par = parammat[4,],
                 fn = function(kappa3){
                   sum( BAMBI::dvmsinmix(data,kappa1 = parammat[2,],
                                         kappa2 = parammat[3,],
                                         kappa3 = kappa3,
                                         mu1 = parammat[5,],
                                         mu2 = parammat[6,],
                                         pmix = parammat[1,],log = TRUE))
                 },
                 lower = -lbound,
                 upper = lbound,
                 method = "L-BFGS-B",
                 control = list(fnscale = -1))
      parammat[4,] <- a$par




    }

    parammat[is.na(parammat)] <- 0.001
    param.seq <- rbind(param.seq,as.vector(parammat))

    complete.data.loglkhd.now <- sum( BAMBI::dvmsinmix(data,kappa1 = parammat[2,],
                                                       kappa2 = parammat[3,],
                                                       kappa3 = parammat[4,],
                                                       mu1 = parammat[5,],
                                                       mu2 = parammat[6,],
                                                       pmix = parammat[1,],log = TRUE))
    loglkhd.seq <- c(loglkhd.seq, complete.data.loglkhd.now)

    # STOP CRITERION ----------------------------------------------------------


    diff <- sum( (param.seq[nrow(param.seq),] - param.seq[nrow(param.seq)-1,])^2, na.rm = TRUE )
    cat(cnt >= maxiter | diff < THRESHOLD)
    if(cnt >= maxiter | diff < THRESHOLD){break}

  }

  list(parammat = parammat, pimat = pimat,
       param.seq = param.seq,
       cnt = cnt,
       loglkhd.seq = loglkhd.seq)


}

sinvM.ECMEb <- function(wstat, kappa1, kappa2, lambda, mu1, mu2, THRESHOLD = 1e-10){

  mu1 <- c(cos(mu1),sin(mu1))
  mu2 <- c(cos(mu2),sin(mu2))

  r <- 0
  while(TRUE){
    mu1c <- kappa1 * wstat$y1bar + lambda * wstat$S12 %*% mu2
    mu1c <- mu1c / sqrt(sum(mu1c^2))
    mu2c <- kappa2 * wstat$y2bar + lambda * t(wstat$S12) %*% mu1
    mu2c <- mu2c / sqrt(sum(mu2c^2))

    # STOP IF NOT CHANGING
    diffb <- norm(mu1c - mu1) + norm(mu2c - mu2)
    r <- r + 1
    if(r >= 100 | diffb < THRESHOLD){break}
    mu1 <- mu1c
    mu2 <- mu2c

  }
  mu1 <- atan2(mu1[2], mu1[1])
  mu2 <- atan2(mu2[2], mu2[1])
  return(list(mu1 = ifelse(mu1 < 0, mu1 + 2*pi, mu1),
              mu2 = ifelse(mu2 < 0, mu2 + 2*pi, mu2)))
}
