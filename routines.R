# final repo for r functions

library(BAMBI)
library(igraph)
library(polynom)

# Dealing with torus data -------------------------------------------------

# on.torus() : transforms data to be on [0, 2pi)^2
# ang.dist() : computes angular distances
# ang.minus(): computes angular subtraction
# ang.pdist(): computes pairwise distances matrix
# grid.torus(): returns an equally-spaced grid on torus
# kde.torus(): returns a kde using circular von Mises

on.torus <- function(x){
  y <- x
  y[,1] <- x[,1] - 2*pi*floor(x[,1]/2/pi)
  y[,2] <- x[,2] - 2*pi*floor(x[,2]/2/pi)
  return(y)
}
ang.dist <- function(x,y){
  # dist(x,y) if x and y are [0, 2*pi]
  apply((rbind( abs(x - y), 
                x + 2*pi - y,
                y + 2*pi - x)),2,min)
} 
ang.minus <- function(x,y){
  # x "minus" y if x and y are [0, 2*pi]
  t <- rbind(x-y, x-y+2*pi, x-y-2*pi)
  
  tind <- apply( abs(t),2, which.min)
  
  tt <- t[1,]
  tt[tind == 2] <- t[2,tind == 2]
  tt[tind == 3] <- t[3,tind == 3]
  
  tt
} 
tor.minus <- function(data, mu){
  # data is n x 2 matrix of toroidal values, mu is a 2-vector. 
  # returns toroidal subtraction 
  cbind(
    ang.minus(data[,1], mu[1]),
    ang.minus(data[,2], mu[2]))
}

ang.pdist <- function(data){
  # assuming that data are n x d angular data on [0,2pi]^d
  # computes L2 angular distance 
  
  n <- nrow(data)
  d <- ncol(data)
  pdistmat <- matrix(0,ncol = n, nrow = n)
  
  for (i in 1:n-1){
    ad <- 0 
    for (dd in 1:d){
      ad <- ad + ( ang.dist(data[i,dd],data[(i+1):n,dd]) )^2
    }
    ad <- sqrt(ad)
    pdistmat[i,(i+1):n] <- ad
    pdistmat[(i+1):n,i] <- ad
  }
  as.dist(pdistmat)
}
grid.torus <- function(grid.size = 100){
  # returns grid points on torus of size (grid.size) x (grid.size)
  Axis <- seq(0, 2*pi, length = grid.size)
  return( cbind(rep(Axis, grid.size),rep(Axis, each=grid.size)))
}
kde.torus <- function(data, eval.point = grid.torus(), 
                      concentration = 25){
  # Computes kde using independent bivariate von mises kernel.
  # evaluated at eval.point, a N x 2 matrix of points on torus
  # returns N-vector of kdes evaluated at eval.point
  # used in cp.torus.kde()
  
  n <- nrow(data)
  N <- nrow(eval.point)
  
  summand <- rep(0,N) 
  for (i in 1:n){
    summand <- summand + exp(concentration * 
                               rowSums( cos( sweep(eval.point,2,data[i,],"-")) ))  
  }
  return( summand / n / (2 * pi * besselI(concentration,0))^2 ) 
}



# conformal prediction and related functions ------------------------------
cp.torus.kde <- function(data, eval.point = grid.torus(), 
                         level= 0.1,
                         concentration = 25){
  # Computes conformal prediction set indices (TRUE if in the set) using kde as conformal scores
  # level can be either a scalar or a vector, or even null.
  # if level is null, return kde at eval.point and at data points.
  # if level is a vector, return the above and Prediction set indices for each value of level. 
  
  N <- nrow(eval.point)
  n <- nrow(data)
  
  eval.point.bind <- rbind(eval.point,data) 
  # 
  phat <- kde.torus(data, eval.point.bind, 
                    concentration = concentration)
  
  phat.grid <- phat[1:N]                 # hat{p}(eval.point)     
  phat.data <- sort( phat[(N+1):(N+n)], index.return = TRUE) 
  data <- data[phat.data$ix,] # data reordered to satisfy y_i = y_(i)
  phat.data <- phat.data$x    # hat{p}(y_(i)) sorted (increasing)
  
  nalpha <- length(level) 
  cp.torus <- NULL
  
  if (!is.null(level)){
    for (i in 1:nalpha){
      ialpha <- floor( (n + 1) * level[i])
      # indices for inclusion in L-
      Lminus <- phat.grid >= phat.data[ialpha]
      
      # indices for inclusion in L+
      const_vM2 <- (2*pi*besselI(concentration,0))^2
      zeta <- (exp(concentration*2) - exp(-concentration*2)) / const_vM2
      Lplus <- phat.grid >=  phat.data[ialpha] - zeta / n   
      
      # indices for inclusion in Cn(alpha)
      K <- (exp(concentration*2) - 
              exp(concentration*(  rowSums(cos(sweep(eval.point,2,data[ialpha,],"-"))  )))) / 
        const_vM2 # K(0) - K(y_(ialpha) - y)
      Cn <- phat.grid - phat.data[ialpha] + K/n >= 0 
      
      cp.torus.i <- data.frame(phi = eval.point[,1], psi = eval.point[,2],
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
 

icp.torus.score <- function(data, split.id = NULL,
                            method = c("all","kde","mixture"),
                            mixturefitmethod = c("circular","axis-aligned","general","Bayesian"),
                            param = list(J = 4,concentration = 25)){
  # returns an icp.torus object, containing all values to compute the conformity score.
  
  # Use sample splitting to produce (inductive) conformal prediction sets
  # if split.id is not supplied, then use a random set of size floor(n/2) for phat estimation
  # split.id is a n-vector consiting of values 1 (estimation) and 2 (evaluation)
  
  # param contains the number of components for mixture fitting
  #       and the concentration parameter.
  
  method <- match.arg(method)
  mixfitmethod <- match.arg(mixturefitmethod)
  
  
  # sample spliting; preparing data
  n <- nrow(data)
  if (is.null(split.id)){
    split.id <- rep(2,n)
    split.id[ sample(n,floor(n/2)) ] <- 1
  }
  
  X1 <- data[split.id == 1,]
  X2 <- data[split.id == 2,]
  n2 <- nrow(X2)
  
  # Prepare output 
  icp.torus <- list(kde = NULL, mixture = NULL,n2 = n2, split.id = split.id)
  
  # For each method, use X1 to estimate phat, then use X2 to provide ranks.
  
  
  # 1. kde 
  if (method != "mixture"){ 
    phat <- kde.torus(X1, X2, concentration = param$concentration)
    # phat.X2.sorted <- sort(phat) 
    
    icp.torus$kde$concentration <- param$concentration
    icp.torus$kde$score <- sort(phat)
    icp.torus$kde$X1 <- X1
  }
  
  # 2. mixture fitting
  if (method != "kde"){
    
    icp.torus$mixture$fittingmethod <- mixturefitmethod
    
    if (mixturefitmethod != "Bayesian"){
      vm2mixfit <- EMsinvMmix(X1, J = param$J, parammat = EMsinvMmix.init(data, param$J), 
                              THRESHOLD = 1e-10,
                              type = mixturefitmethod,
                              kmax = 500,
                              verbose = TRUE)
      icp.torus$mixture$fit <- vm2mixfit
    } else {
      #vm2mixfit ## USE BAMBI  
      stop("Bayesian not yet implemented")
    }
    
    # compute phat(X2)
    phat <- BAMBI::dvmsinmix(X2,kappa1 = vm2mixfit$parammat[2,], 
                             kappa2 = vm2mixfit$parammat[3,],
                             kappa3 = vm2mixfit$parammat[4,],
                             mu1 = vm2mixfit$parammat[5,],
                             mu2 = vm2mixfit$parammat[6,],
                             pmix = vm2mixfit$parammat[1,],log = FALSE)
    icp.torus$mixture$score <- sort(phat)
    
    
    # compute phat_max(X2)
    phatj <- phat.eval(X2, vm2mixfit$parammat)
    
    # phatj <- matrix(0,nrow = n2,ncol = param$J)
    # for(j in 1:param$J){
    #   phatj[,j] <- BAMBI::dvmsin(X2, kappa1 = vm2mixfit$parammat[2,j], 
    #                              kappa2 = vm2mixfit$parammat[3,j],
    #                              kappa3 = vm2mixfit$parammat[4,j],
    #                              mu1 = vm2mixfit$parammat[5,j],
    #                              mu2 = vm2mixfit$parammat[6,j],log = FALSE
    #   ) * vm2mixfit$parammat[1,j]
    # }              
    icp.torus$mixture$score_max <- sort( apply(phatj,1,max)  )
    
    # compute parameters for ellipses and phat_e(X2)
    ellipse.param <- norm.appr.param(vm2mixfit$parammat)
    # ellipse.param <- list(mu1 = NULL, mu2=NULL, Sigmainv = NULL,c = NULL)
    # ellipse.param$mu1 <- vm2mixfit$parammat[5,]
    # ellipse.param$mu2 <- vm2mixfit$parammat[6,]
    # for (j in 1:param$J){
    #   kap1 <- vm2mixfit$parammat[2,j]
    #   kap2 <- vm2mixfit$parammat[3,j]
    #   lamb <- vm2mixfit$parammat[4,j]
    #   pi_j <- vm2mixfit$parammat[1,j]
    #   ellipse.param$Sigmainv[[j]] <- 
    #     matrix(c(kap1, rep(lamb,2), kap2), nrow = 2)
    #   ellipse.param$c[j] <- 2*log(pi_j * (kap1*kap2 - lamb^2) )
    # }
    icp.torus$mixture$ellipsefit <- ellipse.param 
    
    ehatj <- ehat.eval(X2, ellipse.param)
    # ehatj <- matrix(0,nrow = n2,ncol = param$J)
    # for(j in 1:param$J){
    #   z <- tor.minus(X2, c(ellipse.param$mu1[j], ellipse.param$mu2[j]) )
    #   S <- ellipse.param$Sigmainv[[j]]
    #   A <- z %*% S 
    #   ehatj[,j] <- -apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]}) + ellipse.param$c[j]
    # }              
    icp.torus$mixture$score_ellipse <- sort( apply(ehatj,1,max) ) 
    
  }
  
  return(icp.torus)
  
}

phat.eval <- function(X, parammat){
  # evaluate pi_j f_j(x). Returns nrow(X) times ncol(parammat) (n x J) matrix. 
  # used for mixture models
  n2 <- nrow(X)
  J <- ncol(parammat)
  phatj <- matrix(0,nrow = n2,ncol = J)
  for(j in 1:J){
    phatj[,j] <- BAMBI::dvmsin(X, kappa1 = parammat[2,j], 
                               kappa2 = parammat[3,j],
                               kappa3 = parammat[4,j],
                               mu1 = parammat[5,j],
                               mu2 = parammat[6,j],log = FALSE
    ) * parammat[1,j]
  }    
  phatj <- ifelse( is.nan(phatj), 0, phatj)
  return(phatj)
}

norm.appr.param <- function(parammat){
  # parameters for elliptical approximation of bivariate von mises
  J <- ncol(parammat)
  
  ellipse.param <- list(mu1 = NULL, mu2=NULL, Sigmainv = NULL,c = NULL)
  ellipse.param$mu1 <- parammat[5,]
  ellipse.param$mu2 <- parammat[6,]
  for (j in 1:J){
    kap1 <- parammat[2,j]
    kap2 <- parammat[3,j]
    lamb <- parammat[4,j]
    pi_j <- parammat[1,j]
    
    # # small adjustment for low concentration 
    # if (kap1 < 10){
    #   lamb <- lamb * (kap1/10)^(1/10)
    #   kap1 <- kap1 * (kap1/10)^(1/10)
    # }
    # if (kap2 < 10){
    #   lamb <- lamb * (kap2/10)^(1/10)
    #   kap2 <- kap1 * (kap2/10)^(1/10)
    # }
    #  
    
    ellipse.param$Sigmainv[[j]] <- 
      matrix(c(kap1, -rep(lamb,2), kap2), nrow = 2)
    ellipse.param$c[j] <- 2*log(pi_j * (kap1*kap2 - lamb^2) )
  }
  
  return(ellipse.param)
}

ehat.eval <- function(X, ellipse.param){
  
  # evaluate ehat_j(x). Returns nrow(X) times ncol(parammat) (n x J) matrix. 
  # ehat(x) is the columnwise maximum of ehat_j(x): apply(ehatj,1,max)
  # used for elliptically-approximated mixture models
  n2 <- nrow(X)
  J <- length(ellipse.param$mu1)
  
  ehatj <- matrix(0,nrow = n2,ncol = J)  
  for(j in 1:J){
    z <- tor.minus(X, c(ellipse.param$mu1[j], ellipse.param$mu2[j]) )
    S <- ellipse.param$Sigmainv[[j]]
    A <- z %*% S 
    ehatj[,j] <- -apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]}) + ellipse.param$c[j]
  }
  ehatj
}

mah.dist2.eval <- function(X, ellipse.param){
  
  # evaluate mahalanobis_distance from x to each ellipses. Returns nrow(X) times ncol(parammat) (n x J) matrix. 
  
  n2 <- nrow(X)
  J <- length(ellipse.param$mu1)
  
  ehatj <- matrix(0,nrow = n2,ncol = J)  
  for(j in 1:J){
    z <- tor.minus(X, c(ellipse.param$mu1[j], ellipse.param$mu2[j]) )
    S <- ellipse.param$Sigmainv[[j]]
    A <- z %*% S 
    ehatj[,j] <-  apply(cbind(A,z), 1, function(a){a[1]*a[3]+a[2]*a[4]})
  }
  ehatj
}


icp.torus.eval <- function(icp.torus, level = 0.1, eval.point = grid.torus()){
  # evaluates Chat_kde, Chat_mix, Chat_max, Chat_e. 
  N <- nrow(eval.point)
  
  n2 <- icp.torus$n2
  nalpha <- length(level) 
  cp <- list(Chat_kde = NULL, Chat_mix = NULL, Chat_max = NULL, Chat_e = NULL,
             level = level,
             phi = eval.point[,1], 
             psi = eval.point[,2])
  
  
  
  if( !is.null( icp.torus$kde ) ){
    
    Chat_kde <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_kde) <- level
    
    phat.grid <- kde.torus(icp.torus$kde$X1, eval.point, concentration = icp.torus$kde$concentration)
    for (i in 1:nalpha){
      ialpha <- floor( (n2 + 1) * level[i])
      
      # indices for inclusion in Chat_kde 
      Chat_kde[,i] <- phat.grid >= icp.torus$kde$score[ialpha]
    }
    cp$Chat_kde <- Chat_kde
  }
  
  if( !is.null(icp.torus$mixture)){
    Chat_mix <- matrix(0, nrow = N, ncol = nalpha)
    colnames(Chat_mix) <- level
    
    Chat_max <- Chat_mix
    Chat_e <- Chat_mix 
    
    phatj <- phat.eval(eval.point, icp.torus$mixture$fit$parammat) 
    ehatj <- ehat.eval(eval.point, icp.torus$mixture$ellipsefit)
    phat_mix <- rowSums(phatj)
    phat_max <- apply(phatj,1,max)
    ehat <- apply(ehatj,1,max)
    
    for (i in 1:nalpha){
      ialpha <- floor( (n2 + 1) * level[i])
      
      Chat_mix[,i] <- phat_mix >= icp.torus$mixture$score[ialpha]
      Chat_max[,i] <- phat_max >= icp.torus$mixture$score_max[ialpha]
      Chat_e[,i]   <-    ehat  >= icp.torus$mixture$score_ellipse[ialpha]
    }
    
    cp$Chat_mix <- Chat_mix
    cp$Chat_max <- Chat_max
    cp$Chat_e <- Chat_e
  }
  
  return(cp)
  
}










# Functions related to univariate and bivariate sine von Mises ------------

# Ainv()   : returns the inverse of A, where A is the ratio of I_1()/ I_0().
# wtd.stat.ang(): a function to compute (weighted) mean and mean resultant length
# wtd.stat.amb(): a function to compute ambient mean and cross-matrix (S12) 
# BAMBI::dvm(): univariate von Mises density 
# BAMBI::dvmsin(): bivariate sine von Mises density 
# BAMBI::dvmsinmix(): bivariate sine von Mises mixture density 

Ainv <- function (x,kmax = 2000){
  k <- ifelse(0 <= x & x < 0.53, 2 * x + x^3 + (5 * x^5)/6, 
         ifelse(x < 0.85, -0.4 + 1.39 * x + 0.43/(1 - x), 
                ifelse(x < 1,1/(x^3 - 4 * x^2 + 3 * x), kmax)))
  ifelse(k > kmax, kmax, k)
}
wtd.stat.ang <- function(data, w){
  # computes weighted extrinsic mean direction and mean resultant length
  # 
  
  # note w is multiplied to each column of the former
  wtd_ext_mean <- colMeans( ( cos(data) + 1i * sin(data) ) * w)
  
  ang <- Arg(wtd_ext_mean)
  ang <- ifelse(ang < 0, ang + 2*pi,ang)
  R <- Mod(wtd_ext_mean)
  
  return(list(Mean = ang, R = R))
}
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
# BAMBI::dvm() 
# BAMBI::dvmsin()
# BAMBI::dvmsinmix()
 




# Clustering functions ----------------------------------------------------

cluster.assign.torus <- function(data, icp.torus, level){
  # clustering by connected components of ellipses
  # 
  # return clustering assignment for data, given icp.torus objects. 
  # clustering assignment is given by 
  # 1) max_j phatj (phatj = pi_j f_j(x)) for max-mixture (not implemented yet)
  # 2) max_k sum_{j in cluster k} phatj for max-mixture  (not implemented yet)
  # 3) max_j ehatj for ellipse-approx
  # 4) max_k sum_{j in cluster k} ehatj for ellipse-approx
  # two options
  # one: every point is assigned to clusters 1:K by either (1)-(4) above
  # two: outliers are assigned to cluster K+1.
  n2 <- icp.torus$n2
  ialpha <- floor( (n2 + 1) * level)
   
  # For max-ellipse ---------------------------------------------------------
  
  t <- icp.torus$mixture$score_ellipse[ialpha]
  cluster.obj <- conn.comp.ellipse(icp.torus$mixture$ellipsefit, t)
  K <- cluster.obj$ncluster
  ehatj <- ehat.eval(data, icp.torus$mixture$ellipsefit)  
  
  ehatj[,cluster.obj$componentid == 0] <- -Inf
  
  maxj.id <- apply(ehatj, 1, which.max)
  cluster.id1 <- cluster.obj$componentid[maxj.id]
  cluster.obj$cluster.id.by.ehat <- cluster.id1
  
  partsum <- matrix(0, nrow = nrow(data), ncol = K)
  for(k in 1:K){
    ifelse(sum( cluster.obj$componentid == k) > 1, 
           partsum[,k] <- rowSums(ehatj[,cluster.obj$componentid == k]),
           partsum[,k] <- ehatj[,cluster.obj$componentid == k])
  }
  cluster.obj$cluster.id.by.partialsum <- apply(partsum, 1, which.max)
  
  ehat <- apply(ehatj,1,max)
  
  cluster.id1[!(ehat >= icp.torus$mixture$score_ellipse[ialpha])] <- K+1
  cluster.obj$cluster.id.outlier <- cluster.id1
  
  # Use mahalanobis distance 
  
  mah <- mah.dist2.eval(data, icp.torus$mixture$ellipsefit)  
  mah[,cluster.obj$componentid == 0] <- Inf
  
  maxj.id <- apply(mah, 1, which.min)
  cluster.id1 <- cluster.obj$componentid[maxj.id]
  cluster.obj$cluster.id.by.Mah.dist <- cluster.id1
  
  
  return(cluster.obj)
}

conn.comp.ellipse <- function(ellipse.param, t){
  
  require(igraph)
  
  J <- length(ellipse.param$mu1)
  Adj.matrix <- matrix(0,nrow = J, ncol = J)
  for (j in 1:(J-1)){
    for (k in (j+1):J){
      overlap <- Test.intersection.ellipse.torus(ellipse.param, index = c(j,k),t = t)
      Adj.matrix[j,k] <- overlap 
      # Adj.matrix[k,j] <- overlap 
    }
  }
  
  emptysetind <- vector(mode = "logical", length = J)
  for (j in 1:J){
    ifelse( ellipse.param$c[j] - t <= 0 , emptysetind[j] <- TRUE, emptysetind[j] <- FALSE) 
  }
  Graph <- igraph::graph_from_adjacency_matrix(Adj.matrix)
  comp <- igraph::components(Graph)$membership
  comp[emptysetind] <- 0
  comp.id <- unique(comp)
  comp.id <- comp.id[comp.id>0]
  
  K <- length(comp.id)
  conncomp.ind <- comp
  for ( k in 1:K ){
    conncomp.ind[comp == comp.id[k]] <- k 
  }
  
  return(list(componentid = conncomp.ind, ncluster = K))
}

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

Test.intersection.ellipse <- function(mean.1, M.1, mean.2, M.2){
  require(polynom)
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


# Mixture model fitting ---------------------------------------------------

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
    r <- r+1
    if(r >= 100 | diffb < THRESHOLD){break}
    mu1 <- mu1c
    mu2 <- mu2c
    
  }
  mu1 <- atan2(mu1[2],mu1[1])
  mu2 <- atan2(mu2[2],mu2[1]) 
  return(list(mu1 = ifelse(mu1 < 0, mu1 + 2*pi,mu1), 
          mu2 = ifelse(mu2 < 0, mu2 + 2*pi,mu2)))
}
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
     # cat(cnt >= maxiter | diff < THRESHOLD)
    if(cnt >= maxiter | diff < THRESHOLD){break}
    
  }
  
  list(parammat = parammat, pimat = pimat, 
       param.seq = param.seq,
       cnt = cnt,
       loglkhd.seq = loglkhd.seq)
  
  
}

EMsinvMmix.init2 <- function(data, J = 4, hc = hclust(ang.pdist(data))){
  # initialize for EMsinvMmix
  # Use complete-linkage hierarchical clustering
  membership <- cutree(hc, J)
  
  # then use the membership to set mu, kappa by group-wise MLEs. 
  # set the association parameter as zero.
  require(circular)
  
  n <- nrow(data)
  pimat <- matrix(0,n,J)
  for (i in 1:n) {pimat[i,membership[i]] <- 1}
  parammat <- matrix(0,nrow = 6, ncol = J)
  rownames(parammat) <- c("pmix","kappa1","kappa2","kappa3","mu1","mu2")
  parammat[1,] <- colMeans(pimat)
  for(j in 1:J) { 
    a1 <- mle.vonmises(circular(data[membership == j,1]))
    a2 <- mle.vonmises(circular(data[membership == j,2]))
    parammat[2:6,j] <- c(a1$kappa, a2$kappa,0,a1$mu+ pi, a2$mu+pi)
  }
  parammat[is.infinite(parammat)] <- 100
  return(parammat)
}

EMsinvMmix.init <- function(data, J = 4, hc = hclust(ang.pdist(data))){
  # initialize for EMsinvMmix
  # Use complete-linkage hierarchical clustering
  membership <- cutree(hc, J)
  
  # then use the membership to set mu, kappa by group-wise MLEs. 
  # set the association parameter as zero.
  # require(circular)
  kmax <- 500
  n <- nrow(data)
  pimat <- matrix(0,n,J)
  for (i in 1:n) {pimat[i,membership[i]] <- 1}
  parammat <- matrix(0,nrow = 6, ncol = J)
  parammat[1,] <- colMeans(pimat)
  
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
   
  
  rownames(parammat) <- c("pmix","kappa1","kappa2","kappa3","mu1","mu2")
   
  return(parammat)
}
