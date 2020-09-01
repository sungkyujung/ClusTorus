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
