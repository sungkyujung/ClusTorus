# Initalize for EMsinvMmix
#
# \code{EMsinvMmix.init} returns the inital parameter for the optimizing
#   algorithm, \code{EMsinvMmix}.
#
# @details This function uses complete-linkage hierarchical clustering
#   method to set up the initial clusters, where the dissimilarity
#   matrix may be constructed with \code{\link[ClusTorus]{ang.pdist}},
#   which is the distance matrix for angular data. Then, use the
#   membership to set mu, kappa by group-wise MLEs.
#
# @param data n x 2 matrix of toroidal data on \eqn{[0, 2\pi)^2}
# @param J number of components of mixture density
# @param hc use \code{\link[stats]{hclust}} to construct
#   hierarchical clustering. You may tune the parameters of
#   \code{\link[stats]{hclust}} to change the method for
#   initalizing the parameters.
# @return returns \code{parammat} which will become the input of
#   \code{\link[ClusTorus]{EMsinvMmix}}, the detail is below:
#   6 x J parameter data with the following components:
#
#   \code{parammat[1, ]} : the weights for each von Mises sine density
#
#   \code{parammat[n + 1, ]} : \eqn{\kappa_n} for each von Mises
#   sine density for n = 1, 2, 3
#
#   \code{parammat[m + 4, ]} : \eqn{\mu_m} for each von Mises
#    sine density for m = 1, 2
# @references 'S. Jung, K. Park, and B. Kim (2020),
#   "Clustering on the torus by conformal prediction"
# @seealso See the detail for hierarchical clustering:
#   \code{\link[stats]{hclust}}, \code{\link[ClusTorus]{ang.pdist}}.
#    Also, see the detail of \code{\link[ClusTorus]{EMsinvMmix}}.
# @examples
# \dontrun{
# ## mean vectors
#
# Mu1 <- c(3, 0)
# Mu2 <- c(2, 2)
# Mu3 <- c(1, 4)
#
# ## covariance matrices
#
# Sigma1 <- matrix(c(0.1, 0.05, 0.05, 0.2), 2, 2)
# Sigma2 <- matrix(c(0.1, 0, 0, 0.01), 2, 2)
# Sigma3 <- matrix(c(0.01, 0, 0, 0.1), 2, 2)
#
# ## 2-dimensional multivariate normal data wrapped with toroidal space
# require(MASS)
# data <- rbind(mvrnorm(n=70, Mu1, Sigma1),
#               mvrnorm(n=50, Mu2, Sigma2),
#               mvrnorm(n=50, Mu3, Sigma3))
# data <- on.torus(data)
#
# EMsinvMmix.init(data, J = 3)
# }

EMsinvMmix.init <- function(data, J = 4, hc = stats::hclust(ang.pdist(data))){
  # initialize for EMsinvMmix
  # Use complete-linkage hierarchical clustering
  membership <- stats::cutree(hc, J)

  # then use the membership to set mu, kappa by group-wise MLEs.
  # set the association parameter as zero.
  # require(circular)
  kmax <- 500
  n <- nrow(data)
  pimat <- matrix(0, n, J)
  for (i in 1:n) {pimat[i, membership[i]] <- 1}
  parammat <- matrix(0, nrow = 6, ncol = J)
  parammat[1, ] <- colMeans(pimat)

  # For each j, compute weighted angluar mean and a joint mean resultant length
  for(j in 1:J){
    pivec <- pimat[, j]
    wstat <-wtd.stat.ang(data, w = pivec)
    kappa_next <- min(c(Ainv(mean(wstat$R) / parammat[1,j])  , kmax))
    parammat[2, j] <- kappa_next
    parammat[3, j] <- kappa_next
    parammat[4, j] <- 0
    parammat[5:6, j] <- wstat$Mean
  }


  rownames(parammat) <- c("pmix", "kappa1", "kappa2", "kappa3", "mu1", "mu2")

  return(parammat)
}

#
# EMsinvMmix.init2 <- function(data, J = 4, hc = hclust(ang.pdist(data))){
#   # initialize for EMsinvMmix
#   # Use complete-linkage hierarchical clustering
#   membership <- cutree(hc, J)
#
#   # then use the membership to set mu, kappa by group-wise MLEs.
#   # set the association parameter as zero.
#   require(circular)
#
#   n <- nrow(data)
#   pimat <- matrix(0,n,J)
#   for (i in 1:n) {pimat[i,membership[i]] <- 1}
#   parammat <- matrix(0,nrow = 6, ncol = J)
#   rownames(parammat) <- c("pmix","kappa1","kappa2","kappa3","mu1","mu2")
#   parammat[1,] <- colMeans(pimat)
#   for(j in 1:J) {
#     a1 <- mle.vonmises(circular(data[membership == j,1]))
#     a2 <- mle.vonmises(circular(data[membership == j,2]))
#     parammat[2:6,j] <- c(a1$kappa, a2$kappa,0,a1$mu+ pi, a2$mu+pi)
#   }
#   parammat[is.infinite(parammat)] <- 100
#   return(parammat)
# }
