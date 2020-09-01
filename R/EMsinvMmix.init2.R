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
