# Connected components of given ellipses
#
# \code{conn.comp.ellipse} investigates how many connected
#   components exist and which ellipses are included in
#   each connected component.
#
# @param ellipse.param list which is consisting of mean of each angular
#   coordinate, inverse of each covariance matrix, and constant term.
#
# @param t a numeric value which determines the size of ellipses.
# @return a \code{list} which contains two components:
#
#   \code{componentid} indicates which connected component
#   the ellipse belongs to.
#   \code{ncluster} indicates the number of clusters, i.e.,
#   the number of connected components.
# @references 'S. Jung, K. Park, and B. Kim (2020),
#   "Clustering on the torus by conformal prediction"
# @examples
# \dontrun{
#
# parammat <- matrix(c(0.4, 0.3, 0.3,
#                      20, 25, 25,
#                      30, 25, 20,
#                      1, 2, 3,
#                      1, 2, 3,
#                      0, 2, 4), nrow = 6, byrow =TRUE)
#
# ellipse.param <- norm.appr.param(parammat)
#
# t <- 0.5
#
# conn.comp.ellipse(ellipse.param, t)
# }

conn.comp.ellipse <- function(ellipse.param, t){

  J <- length(ellipse.param$c)
  d <- ncol(ellipse.param$Sigmainv[[1]])
  Adj.matrix <- matrix(0,nrow = J, ncol = J)

  combinations <- t(utils::combn(1:J, 2))
  # for (j in 1:(J-1)){
  #   for (k in (j+1):J){
  #     overlap <- Test.intersection.ellipse.torus(ellipse.param, index = c(j,k),t = t)
  #     Adj.matrix[j,k] <- overlap
  #     # Adj.matrix[k,j] <- overlap
  #   }
  # }
  overlap.results <- purrr::map_int(1:nrow(combinations),
                                    function(i) {Test.intersection.ellipse.torus(ellipse.param, index = combinations[i, ],t = t)})

  for (i in 1:nrow(combinations)){
    Adj.matrix[combinations[i, 1], combinations[i, 2]] <- overlap.results[i]
  }

  emptysetind <- vector(mode = "logical", length = J)
  for (j in 1:J){
    ifelse( (ellipse.param$c[j] - t <= 0) || (det(ellipse.param$Sigmainv[[j]]) >= 1e+6^d) ,
            emptysetind[j] <- TRUE, emptysetind[j] <- FALSE)
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
