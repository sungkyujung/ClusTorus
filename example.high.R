library(MASS)
library(tidyverse)
library(ClusterR)
library(mclust)
library(cowplot)
# library(ClusTorus)
require(GGally)
devtools::load_all()

set.seed(20201)
# Prepare data 1 ------------------------------------------------------------
Mu1 <- c(3,0,0,0,0,0,0)
Mu2 <- c(2,2,2,2,2,2,2)
Mu3 <- c(1,4,4,4,4,4,4)
# Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
Sigma1 <- matrix(rep(0.05, 49), nrow = 7) - diag(0.05, 7) + diag(c(0.1,0.2,0.1,0.2,0.1,0.2,0.1))
# Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
Sigma2 <- diag(c(0.1,0.01,0.1,0.01,0.1,0.01,0.1))
# Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)
Sigma3 <- diag(c(0.01,0.1,0.01,0.1,0.01,0.1,0.01))

Example1 <- rbind(mvrnorm(n=70, Mu1, Sigma1),
                  mvrnorm(n=50, Mu2, Sigma2),
                  mvrnorm(n=50, Mu3, Sigma3)
                  # data.unif,
                  # data.diamond
                  )

Example1 <- on.torus(Example1)

label <- c(rep(1,70),rep(2,50),
           rep(3,50)
           # rep(4,50),
           # rep(5,50)
           )

dat1 <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat1) <- c("phi","psi","chi1","chi2","chi3","chi4","chi5","label")

Example1 <- rbind(mvrnorm(n=70, Mu1, Sigma1),
                  mvrnorm(n=50, Mu2, Sigma2),
                  mvrnorm(n=50, Mu3, Sigma3)
                  # data.unif,
                  # data.diamond
)

Example1 <- on.torus(Example1)

label <- c(rep(1,70),rep(2,50),
           rep(3,50)
           # rep(4,50),
           # rep(5,50)
)

dat1.test <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat1.test) <- c("phi","psi","chi1","chi2","chi3","chi4","chi5","label")

# Clustering by four methods (subfunction definition)------------------------------------------
Example_paper_supp <- function(J, dat1, dat1.test, type = c("homogeneous-circular",
                                                            "heterogeneous-circular",
                                                            "ellipsoids",
                                                            "general")){

  type <- match.arg(type)
  data <- dat1[,1:7]
  data.test <- dat1.test

  data.test.label <- data.test[,8]
  data.test <- data.test[,1:7]

  predicted.label <- matrix(NA, nrow = nrow(data), ncol =2)
  # colnames(predicted.label) <- c("kmeans_naive","kmeans_amb","Proposal_ehat","Proposal_out")
  #
  # # Kmeans_naive
  # kmeans.out<-KMeans_rcpp(data,clusters = J)
  # predicted.label[,1] <- predict_KMeans(data.test, kmeans.out$centroids)
  #
  # # kmeans_ambient space
  # kmeans.out<-KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
  # predicted.label[,2] <- predict_KMeans(cbind(cos(data.test),sin(data.test)) , kmeans.out$centroids)

  # Clustering by our method
  # start <- Sys.time()
  # 1) Find alpha and J
  Jvec <- 3:35
  l <- list()

  # sample spliting; preparing data
  set.seed(2020)
  n <- nrow(data)
  split.id <- rep(2,n)
  split.id[ sample(n,floor(n/2)) ] <- 1

  for (j in Jvec){
    l[[j]] <- icp.torus.score(as.matrix(data), split.id = split.id,
                              method = "kmeans",
                              kmeansfitmethod = type,
                              init = "k",
                              additional.condition = T,
                              param = list(J = j))
  }

  n2 <- l[[10]]$n2
  alphavec <- 1:floor(n2/2) / n2
  N <- length(alphavec)
  out <- data.frame()

  # end <- Sys.time()
  # cat (end - start)
  # cat("\n")

  # start <- Sys.time()
  # for (j in Jvec){
  #   Mvec <- alphavec
  #   a<-icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus(d = 7, grid.size = 10))
  #   # for (i in 1:N){
  #   #   Mvec[i] <- sum(a$Chat_kmeans[, i])/10000
  #   # }
  #   Mvec <- colSums(a$Chat_kmeans)/10^7
  #
  #   out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec +  Mvec))
  #
  # }
  Jhat <- 10
  alphahat <- 0.05
  out <- rbind(out, data.frame(alpha = alphahat, J = Jhat))

  end <- Sys.time()
  cat (end - start)
  cat("\n")

  out.index <- which.min(out$criterion)
  out[out.index,]

  Jhat <- out[out.index,2]
  alphahat <- out[out.index,1]
  icp.torus <- l[[Jhat]]

  ia <- icp.torus.eval(icp.torus, level = alphahat, eval.point = grid.torus(d = 7, grid.size = 10))
  b <- data.frame(ia$eval.point, ia$Chat_kmeans == 1)
  colnames(b) <- c("phi","psi", "C_kmeans")
  head(b)

  # b<- b %>%  pivot_longer(3, names_to = "Type", values_to = "value")
  # g0 <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b, size = 1,lineend = "round" ) +
  #   geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) +
  #   scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
  # g0 <- g0 + ggtitle(paste("ICP, kellipsoids fitting with J=",Jhat,", alpha=",alphahat))
  # g0

  #c <- cluster.assign.torus(data, icp.torus, level = alphahat)
  #c
  start <- Sys.time()
  c <- cluster.assign.torus(data.test, icp.torus, level = alphahat)
  end <- Sys.time()

  cat (end - start)
  cat("\n")
  c
  predicted.label[,1] <- c$kmeans$cluster.id.by.ehat
  predicted.label[,2] <- c$kmeans$cluster.id.outlier

  aa <- rep(0,2)

  for (j in 1:2){
    aa[j] <- adjustedRandIndex(predicted.label[,j],data.test.label)
  }

  # g_e <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(c$kmeans$cluster.id.by.ehat)) %>%
  #   ggplot(aes(phi,psi, color = membership)) + geom_point() +
  #   scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   ggtitle(paste("Clustering by kmeans to kspheres, K=",c$kmeans$ncluster))
  #
  #
  # g1 <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(c$kmeans$cluster.id.outlier)) %>%
  #   mutate(membership = ifelse(membership == max(c$kmeans$cluster.id.outlier), "out", membership)) %>%
  #   ggplot(aes(phi,psi, color = membership)) + geom_point() +
  #   scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   ggtitle(paste("Clusters and outliers, K=",length(unique(c$kmeans$cluster.id.outlier))))
  #
  #
  # g_kmeans1 <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(predicted.label[,1])) %>%
  #   ggplot(aes(phi,psi, color = membership)) + geom_point() +
  #   scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   ggtitle(paste("Ordinary Kmeans , K=",J))
  #
  #
  # g_kmeans2 <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(predicted.label[,2])) %>%
  #   ggplot(aes(phi,psi, color = membership)) + geom_point() +
  #   scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  #   ggtitle(paste("Kmeans (ambient), K=",J))


  ######################## Empirical coverage
  #
  # start <- Sys.time()

  grid.test <- grid.torus()
  testing.n <- nrow(data)
  data.test <-as.matrix(data.test)
  # given a  C, for each testing item, see if it is included in C.

  # aggregate C's

  alphavec <- 1:floor(n2/2) / n2
  alphavec <- alphavec[alphavec <= 0.4]
  N <- length(alphavec)
  CC <- icp.torus.eval(icp.torus, level = alphavec)

  C <- CC$Chat_kmeans
  Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  for (j in 1:testing.n){
    Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  }
  Out.dt <- data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_kmeans")

  #   icp.torus.kde<- icp.torus.score(as.matrix(data), split.id = split.id,
  #                                   method = "kde",
  #                                   mixturefitmethod = "a",
  #                                   param = list(concentration = 25))
  # #
  #   icp.kde.region <- icp.torus.eval(icp.torus.kde, level = alphavec, eval.point = grid.torus())
  #
  #   C <- icp.kde.region$Chat_kde
  #   Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  #   for (j in 1:testing.n){
  #     Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  #   }
  #   Out.dt <- rbind(Out.dt,
  #                   data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_KDE"))
  #

  head(Out.dt)
  g_cover <- Out.dt %>% ggplot(aes(1-alpha, coverage,color = Type, group = Type))  + geom_line() +
    geom_point() + geom_abline(slope = 1, intercept = 0) + geom_vline(xintercept =  1- alphahat,linetype = "dotted" )



  # end <- Sys.time()
  # cat (end - start)
  # cat("\n")

  ##########################







  list(choice = out[out.index,], labels = predicted.label, Rand = aa, clustering.result = c,
       gg = g0,ge = g_e, gout = g1, g_k1 = g_kmeans1, g_k2 = g_kmeans2, coverage = Out.dt, g_cover  = g_cover )

}

# prespecified J and alpha ---------------------------------------------
# type <- match.arg(type)
type <- "he"

# predetermine J, alpha
Jhat <- 10
alphahat <- 0.05

data <- dat1[, 1:7]
data.test <- dat1.test

data.test.label <- data.test[, 8]
data.test <- data.test[, 1:7]

predicted.label <- matrix(NA, nrow = nrow(data), ncol = 2)

Jvec <- 3:35
l <- list()

# sample spliting; preparing data
set.seed(2020)
n <- nrow(data)
split.id <- rep(2,n)
split.id[ sample(n,floor(n/2)) ] <- 1

# testing for each J
for (j in Jvec){
  l[[j]] <- icp.torus.score(as.matrix(data), split.id = split.id,
                            method = "kmeans",
                            kmeansfitmethod = type,
                            init = "k",
                            additional.condition = T,
                            param = list(J = j))
}

# -----------------
# testing for various J, alpha
# n2 <- l[[10]]$n2
# alphavec <- 1:floor(n2/2) / n2
# N <- length(alphavec)
out <- data.frame()

# evaluating the Lebesgue measure
# for (j in Jvec){
#   Mvec <- alphavec
#   a<-icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus(d = 7, grid.size = 10))
#   # for (i in 1:N){
#   #   Mvec[i] <- sum(a$Chat_kmeans[, i])/10000
#   # }
#   Mvec <- colSums(a$Chat_kmeans)/10^7
#
#   out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec +  Mvec))
#
# }
# ------------------

Jhat <- 10
alphahat <- 0.05
out <- rbind(out, data.frame(alpha = alphahat, J = Jhat))

# out.index <- which.min(out$criterion)
# out[out.index,]
#
# Jhat <- out[out.index,2]
# alphahat <- out[out.index,1]
icp.torus <- l[[Jhat]]

# ia <- icp.torus.eval(icp.torus, level = alphahat, eval.point = grid.torus(d = 7, grid.size = 10))
# b <- data.frame(ia$eval.point, ia$Chat_kmeans == 1)
# colnames(b) <- c("phi","psi", "C_kmeans")
# head(b)

c <- cluster.assign.torus(data.test, icp.torus, level = alphahat)

predicted.label[,1] <- c$kmeans$cluster.id.by.ehat
predicted.label[,2] <- c$kmeans$cluster.id.outlier

aa <- rep(0,2)

for (j in 1:2){
  aa[j] <- adjustedRandIndex(predicted.label[,j],data.test.label)
}


# Actual run --------------------------------------------------------------

# J <- 5
# dat1.results.ho <- Example_paper_supp(J, dat1, dat1.test, type = "ho")
# dat1.results.he <- Example_paper_supp(J, dat1, dat1.test, type = "he")
# dat1.results.e <- Example_paper_supp(J, dat1, dat1.test, type = "e")
# dat1.results.ge <- Example_paper_supp(J, dat1, dat1.test, type = "g")
# J <- 2
# dat2.results.ho <- Example_paper_supp(J, dat2, dat2.test, type = "ho")
# dat2.results.he <- Example_paper_supp(J, dat2, dat2.test, type = "he")
# dat2.results.e <- Example_paper_supp(J, dat2, dat2.test, type = "e")
# dat2.results.ge <- Example_paper_supp(J, dat2, dat2.test, type = "g")
aa

# GGally::ggpairs(dat1, aes(color = label))

result.dat <- data.frame(data.test, membership = as.factor(c$kmeans$cluster.id.outlier)) %>%
  mutate(membership = ifelse(membership == max(c$kmeans$cluster.id.outlier), "out", membership))

GGally::ggpairs(result.dat[, 1:7], aes(color = result.dat[, 8]))












