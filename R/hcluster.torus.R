hcluster.torus <- function(data, J = 4, method = "complete", members = NULL){

  # hierarchical clustering for the data
  # this function will not be exported: only use for initialzing method for kmeans.kspheres
  n <- nrow(data)
  d <- ncol(data)
  hc = stats::hclust(ang.pdist(data), method = method, members = members)

  membership <- stats::cutree(hc, J)

  # -------- initializing ---------------
  hcluster <- list(data = data, centers = NULL,
                   membership = membership, totss = NULL, withinss = NULL,
                   betweenss = NULL, size = NULL)

  hcluster$size <- rep(0, J)
  hcluster$centers <- matrix(0, nrow = J, ncol = d)
  hcluster$withinss <- rep(0, J)

  # -------- calculate withinss -----------
  for(j in 1:J){
    hcluster$size[j] <- sum(membership == j)
    nj <- hcluster$size[j]
    hcluster$centers[j, ] <- wtd.stat.ang(data[membership == j, ],
                                          w = rep(1, nj) / nj)$Mean

    # if the size of cluster is 1, withinss is 0
    if (nj == 1) { next }

    # the case for the cluster size larger than 1
    j.mean <- hcluster$centers[j, ]
    hcluster$withinss[j] <- sum(tor.minus(data[hcluster$membership == j, ], j.mean)^2)
  }

  # --------- calculate totss ------------

  tot.mean <- wtd.stat.ang(data, rep(1, n) / n)$Mean
  hcluster$totss <- sum(tor.minus(data, tot.mean)^2)

  # --------- calculate betweenss ----------

  # if there is only one cluster, then center.distmat must be 0
  center.distmat <- ang.pdist(hcluster$centers)
  hcluster$betweenss <- sum(center.distmat^2)

  return(hcluster)
}
