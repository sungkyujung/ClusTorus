# naturally dispersed data (Park Ex 2)
library(MASS)
library(tidyverse)
library(ClusterR)
library(mclust)
library(cowplot)
devtools::load_all()
# library(ClusTorus)
# source('routines.R')

set.seed(20201)
# Prepare data 1 ------------------------------------------------------------

# Example: five clusters
Mu1 <- c(3,0)
Mu2 <- c(2,2)
Mu3 <- c(1,4)
Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)

unidata <- cbind(2*runif(50, -0.5, 0.5), 0.5*runif(50, -0.5, 0.5))
data.unif <- cbind(unidata[,1]+ 0, unidata[,2] + 1)
data.diamond <- t(matrix(c(cos(-pi/4),-sin(-pi/4),
                           sin(-pi/4),cos(-pi/4)),2,2) %*% t(unidata)) +cbind(rep(5, 50), rep(3, 50))

Example1 <- rbind(mvrnorm(n=70, Mu1, Sigma1),
                  mvrnorm(n=50, Mu2, Sigma2),
                  mvrnorm(n=50, Mu3, Sigma3),
                  data.unif,
                  data.diamond)

Example1 <- on.torus(Example1)
label <- c(rep(1,70),rep(2,50),
           rep(3,50),
           rep(4,50),
           rep(5,50))
dat1 <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat1) <- c("phi","psi","label")

unidata <- cbind(2*runif(50, -0.5, 0.5), 0.5*runif(50, -0.5, 0.5))
data.unif <- cbind(unidata[,1]+ 0, unidata[,2] + 1)
data.diamond <- t(matrix(c(cos(-pi/4),-sin(-pi/4),
                           sin(-pi/4),cos(-pi/4)),2,2) %*% t(unidata)) +cbind(rep(5, 50), rep(3, 50))

Example1 <- rbind(mvrnorm(n=70, Mu1, Sigma1),
                  mvrnorm(n=50, Mu2, Sigma2),
                  mvrnorm(n=50, Mu3, Sigma3),
                  data.unif,
                  data.diamond)
Example1 <- on.torus(Example1)
label <- c(rep(1,70),rep(2,50),
           rep(3,50),
           rep(4,50),
           rep(5,50))
dat1.test <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat1.test) <- c("phi","psi","label")


# Prepare data 2 ------------------------------------------------------------

set.seed(2020)
# Example: ball and L shape
Mu <- c(1,5)
Sigma <- matrix(c(0.05,0,0,0.05),2,2)

unidata.unif <- cbind(runif(50, 0, 2*pi), runif(50, 0, 2*pi))
data.unif1 <- cbind(1.25*runif(150, -1, 1)+ 3.25, 0.5*runif(150, -1, 1) + 1)
data.unif2 <- cbind(0.5*runif(200, -1, 1)+ 5, 2*runif(200, -1, 1) + 2.5)

Example1 <- rbind(mvrnorm(n=100, Mu, Sigma),
                  unidata.unif,
                  data.unif1,
                  data.unif2)
Example1 <- Example1 + 2

Example1 <- on.torus(Example1)
label <- c(rep(1,100),
           rep(3,50),
           rep(2,350))
dat2 <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat2) <- c("phi","psi","label")
unidata.unif <- cbind(runif(50, 0, 2*pi), runif(50, 0, 2*pi))
data.unif1 <- cbind(1.25*runif(150, -1, 1)+ 3.25, 0.5*runif(150, -1, 1) + 1)
data.unif2 <- cbind(0.5*runif(200, -1, 1)+ 5, 2*runif(200, -1, 1) + 2.5)

Example1 <- rbind(mvrnorm(n=100, Mu, Sigma),
                  unidata.unif,
                  data.unif1,
                  data.unif2)
Example1 <- Example1 + 2

Example1 <- on.torus(Example1)
label <- c(rep(1,100),
           rep(3,50),
           rep(2,350))
dat2.test <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat2.test) <- c("phi","psi","label")


# plot(dat2, xlim=c(0, 2*pi), ylim=c(0, 2*pi), xlab="phi", ylab="psi")


# Figure ------------------------------------------------------------------
g1 <- dat1 %>% ggplot(aes(x = phi, y = psi, color = label)) + geom_point() +
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  ggtitle('Data set 1 with true labels')
g2 <- dat2 %>% ggplot(aes(x = phi, y = psi, color = label)) + geom_point() +
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  ggtitle('Data set 2 with true labels')

plot_grid(g1,g2, label_size = 12)

ggsave("./examples/Toy_Data0.png", width = 8, height = 4)



# Clustering by four methods (subfunction definition)------------------------------------------
Example_paper_supp <- function(J, dat1, dat1.test, type = c("homogeneous-circular",
                                                            "heterogeneous-circular",
                                                            "ellipsoids",
                                                            "general")){

  type <- match.arg(type)
  data <- dat1[,1:2]
  data.test <- dat1.test

  data.test.label <- data.test[,3]
  data.test <- data.test[,1:2]

  predicted.label <- matrix(NA, nrow = nrow(data), ncol =4)
  colnames(predicted.label) <- c("kmeans_naive","kmeans_amb","Proposal_ehat","Proposal_out")

  # Kmeans_naive
  kmeans.out<-KMeans_rcpp(data,clusters = J)
  predicted.label[,1] <- predict_KMeans(data.test, kmeans.out$centroids)

  # kmeans_ambient space
  kmeans.out<-KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
  predicted.label[,2] <- predict_KMeans(cbind(cos(data.test),sin(data.test)) , kmeans.out$centroids)

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

  for (j in Jvec){
    Mvec <- alphavec
    a<-icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus())
    # for (i in 1:N){
    #   Mvec[i] <- sum(a$Chat_kmeans[, i])/10000
    # }
    Mvec <- colSums(a$Chat_kmeans)/10000

    out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec +  Mvec))

  }

  out.index <- which.min(out$criterion)
  out[out.index,]

  Jhat <- out[out.index,2]
  alphahat <- out[out.index,1]
  icp.torus <- l[[Jhat]]

  ia <- icp.torus.eval(icp.torus, level = alphahat, eval.point = grid.torus(d = 2, grid.size = 100))
  b <- data.frame(ia$eval.point, ia$Chat_kmeans == 1)
  colnames(b) <- c("phi","psi", "C_kmeans")
  head(b)

  b<- b %>%  pivot_longer(3, names_to = "Type", values_to = "value")
  g0 <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b, size = 1,lineend = "round" ) +
    geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) +
    scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
  g0 <- g0 + ggtitle(paste("ICP, kellipsoids fitting with J=",Jhat,", alpha=",alphahat))
  g0

  #c <- cluster.assign.torus(data, icp.torus, level = alphahat)
  #c
  # start <- Sys.time()
  c <- cluster.assign.torus(data.test, icp.torus, level = alphahat, intersection.plot = TRUE,
                            coord = t(combn(1:2, 2)))
  # end <- Sys.time()

  # cat (end - start)
  # cat("\n")
  #c
  predicted.label[,3] <- c$kmeans$cluster.id.by.ehat
  predicted.label[,4] <- c$kmeans$cluster.id.outlier

  aa <- rep(0,4)

  for (j in 1:4){
    aa[j] <- adjustedRandIndex(predicted.label[,j],data.test.label)
  }

  g_e <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(c$kmeans$cluster.id.by.ehat)) %>%
    ggplot(aes(phi,psi, color = membership)) + geom_point() +
    scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    ggtitle(paste("Clustering by kmeans to kspheres, K=",c$kmeans$ncluster))


  g1 <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(c$kmeans$cluster.id.outlier)) %>%
    mutate(membership = ifelse(membership == max(c$kmeans$cluster.id.outlier), "out", membership)) %>%
    ggplot(aes(phi,psi, color = membership)) + geom_point() +
    scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    ggtitle(paste("Clusters and outliers, K=",length(unique(c$kmeans$cluster.id.outlier))))


  g_kmeans1 <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(predicted.label[,1])) %>%
    ggplot(aes(phi,psi, color = membership)) + geom_point() +
    scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    ggtitle(paste("Ordinary Kmeans , K=",J))


  g_kmeans2 <- data.frame(phi = data.test[,1], psi = data.test[,2], membership = as.factor(predicted.label[,2])) %>%
    ggplot(aes(phi,psi, color = membership)) + geom_point() +
    scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
    ggtitle(paste("Kmeans (ambient), K=",J))


  ######################## Empirical coverage
  #
  # start <- Sys.time()

  grid.test <- grid.torus(d = 2, grid.size = 100)
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


# Actual run --------------------------------------------------------------

J <- 5
dat1.results.ho <- Example_paper_supp(J, dat1, dat1.test, type = "ho")
dat1.results.he <- Example_paper_supp(J, dat1, dat1.test, type = "he")
dat1.results.e <- Example_paper_supp(J, dat1, dat1.test, type = "e")
dat1.results.ge <- Example_paper_supp(J, dat1, dat1.test, type = "g")
J <- 2
dat2.results.ho <- Example_paper_supp(J, dat2, dat2.test, type = "ho")
dat2.results.he <- Example_paper_supp(J, dat2, dat2.test, type = "he")
dat2.results.e <- Example_paper_supp(J, dat2, dat2.test, type = "e")
dat2.results.ge <- Example_paper_supp(J, dat2, dat2.test, type = "g")

model1 <- cbind(dat1.results.ho$Rand[4], dat1.results.he$Rand[4],
                dat1.results.e$Rand[4], dat1.results.ge$Rand[4])
model2 <- cbind(dat2.results.ho$Rand[4], dat2.results.he$Rand[4],
                dat2.results.e$Rand[4], dat2.results.ge$Rand[4])

model1
model2
# dat1.results$Rand
# dat2.results$Rand

plot_grid(dat1.results.ho$gg ,
          dat2.results.ho$gg ,
          dat1.results.he$gg ,
          dat2.results.he$gg ,
          dat1.results.e$gg ,
          dat2.results.e$gg ,
          dat1.results.ge$gg ,
          dat2.results.ge$gg ,
          label_size = 12, nrow = 4, ncol = 2
          )
ggsave("./examples/Toy_Data1.png", width = 8, height = 4)

plot_grid(
  dat1.results.ho$g_k1 ,
  dat2.results.ho$g_k1 ,
  dat1.results.ho$g_k2 ,
  dat2.results.ho$g_k2 , label_size = 12, nrow = 2, ncol = 2)

plot_grid(
  dat1.results.ho$gout ,
  dat2.results.ho$gout ,
  dat1.results.he$gout ,
  dat2.results.he$gout ,
  dat1.results.e$gout ,
  dat2.results.e$gout ,
  dat1.results.ge$gout ,
  dat2.results.ge$gout ,
  label_size = 12, nrow = 4, ncol = 2)

ggsave("./examples/Toy_Data2.png", width = 12, height = 16)

# dat1.results$Rand
# dat2.results$Rand

l1<- dat1.results$coverage %>% ggplot(aes(1-alpha, coverage,color = Type, group = Type))  + geom_line() +
  geom_point() + geom_abline(slope = 1, intercept = 0)
l2 <-dat2.results$coverage %>% ggplot(aes(1-alpha, coverage,color = Type, group = Type))  + geom_line() +
  geom_point() + geom_abline(slope = 1, intercept = 0)


plot_grid(l1 ,
          l2 , label_size = 12)
ggsave("./examples/Toy_Data3.png", width = 8, height = 4)

