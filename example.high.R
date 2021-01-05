library(MASS)
library(tidyverse)
library(ClusterR)
library(mclust)
library(cowplot)
# library(ClusTorus)
require(GGally)
devtools::load_all()

# Prepare data 1 (4-dimensional case) ------------------------------------------------------------
Mu1 <- c(0.89, -0.34, 1.08, 3.02)
Mu2 <- c(1.94, 2.39, -1.10, -1.07)
Mu3 <- c(1.57, -0.97, -1.00, -1.06)
Mu4 <- c(1.48, -1.00, -1.06, 3.00)
Mu5 <- c(1.04, -0.92, -1.05, 2.94)
Mu6 <- c(2.04, 2.36, -1.17, 2.94)


# Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
Sigma.inv <- matrix(2, nrow = 4, ncol = 4) + diag(8, nrow = 4)
Sigma <- solve(Sigma.inv)
# # Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
# Sigma2 <- diag(c(0.1,0.01,0.1,0.01,0.1,0.01,0.1))
# # Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)
# Sigma3 <- diag(c(0.01,0.1,0.01,0.1,0.01,0.1,0.01))

Example1 <- rbind(mvrnorm(n = 125, Mu1, Sigma),
                  mvrnorm(n = 125, Mu2, Sigma),
                  mvrnorm(n = 125, Mu3, Sigma),
                  mvrnorm(n = 125, Mu4, Sigma),
                  mvrnorm(n = 250, Mu5, Sigma),
                  mvrnorm(n = 250, Mu6, Sigma)
                  # data.unif,
                  # data.diamond
                  )

Example1 <- on.torus(Example1)

label <- c(rep(1,125),rep(2,125),
           rep(3,125),rep(4,125),
           rep(4,250),rep(5,250)
           # rep(4,50),
           # rep(5,50)
           )

dat1 <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat1) <- c("phi","psi","chi1","chi2","label")

Example1 <- rbind(mvrnorm(n = 125, Mu1, Sigma),
                  mvrnorm(n = 125, Mu2, Sigma),
                  mvrnorm(n = 125, Mu3, Sigma),
                  mvrnorm(n = 125, Mu4, Sigma),
                  mvrnorm(n = 250, Mu5, Sigma),
                  mvrnorm(n = 250, Mu6, Sigma)
                  # data.unif,
                  # data.diamond
                  )

Example1 <- on.torus(Example1)

label <- c(rep(1,125),rep(2,125),
           rep(3,125),rep(4,125),
           rep(4,250),rep(5,250)
           # rep(4,50),
           # rep(5,50)
           )

dat1.test <- cbind( as.data.frame(Example1) , as.factor(label))
colnames(dat1.test) <- c("phi","psi","chi1","chi2","label")

GGally::ggpairs(dat1[, 1:4], aes(color = dat1[, 5]))

# prespecified J and alpha ---------------------------------------------
# type <- match.arg(type)
type <- "he"

# predetermine J, alpha

data <- dat1[, -ncol(dat1)]
data.test <- dat1.test

data.test.label <- data.test[, ncol(dat1.test)]
data.test <- data.test[, -ncol(data.test)]

predicted.label <- matrix(NA, nrow = nrow(data), ncol = 2)

Jvec <- 10:80
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

cluster.num <- cluster.assign.number(data, Jmin = 10, Jmax = 80)
cluster.num$plot
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

Jhat <- 20
alphahat <- 0.1
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
# start <- Sys.time()
c <- cluster.assign.torus(data.test, icp.torus, level = alphahat)
# end <- Sys.time()

# cat(end - start)

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

# GGally::ggpairs(result.dat[, 1:7], aes(color = result.dat[, 8]))
K <- length(unique(result.dat[,5]))
label <- as.factor(result.dat[, 5])
pairs(result.dat[, 1:4], pch = 19, col = rainbow(K)[label], upper.panel = NULL)
# ggpairs(result.dat[, 1:4], aes(color = result.dat[, 5]))

# case1 : extrinsic kmeans
J <- 5

kmeans.out<-KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
predicted.label1 <- predict_KMeans(cbind(cos(data.test),sin(data.test)) , kmeans.out$centroids)

adjustedRandIndex(predicted.label1,data.test.label)

# case2 : hierarchical
h.out <- hclust(d = ang.pdist(data), method = "complete")
predicted.label2 <- cutree(h.out, J)

adjustedRandIndex(predicted.label2,data.test.label)
# kmeans.out2 <- kmeans.torus(data, J)
# predicted.label2 <- predict.kmeans.torus(data.test, kmeans.out2)
# adjustedRandIndex(predicted.label2,data.test.label)


set.seed(2021)
# Prepare data 2 (6-dimensional case) ------------------------------------------------------------
Mu1 <- c(-0.94, -1.94, 0.86, 0.86, 0.86, 0.86)
Mu2 <- c(-3.14, -3.14, -3.14, -3.14, -3.14, -3.14)
Mu3 <- c(-0.74, 0.06, 1.36, 1.86, 1.86, 1.86)
Mu4 <- c(-0.2, -0.2, -0.2, -0.2, 0.66, 0.66)
Mu5 <- c(-0.64, -1.64, -0.64, -1.64, -0.64, -1.64)
Mu6 <- c(1.14, 1.14, 1.14, 1.14, 1.14, 1.14)


# Sigma1 <- matrix(c(0.1,0.05,0.05,0.2),2,2)
Sigma1.inv <- matrix(4.5, nrow = 6, ncol = 6) + diag(15.5, nrow = 6)
Sigma2.inv <- matrix(-3, nrow = 6, ncol = 6) + diag(23, nrow = 6)
Sigma3.inv <- matrix(-2, nrow = 6, ncol = 6) + diag(22, nrow = 6)
Sigma4.inv <- matrix(5, nrow = 6, ncol = 6) + diag(15, nrow = 6)
Sigma5.inv <- matrix(2, nrow = 6, ncol = 6) + diag(18, nrow = 6)
Sigma6.inv <- matrix(-2, nrow = 6, ncol = 6) + diag(22, nrow = 6)
Sigma1 <- solve(Sigma1.inv)
Sigma2 <- solve(Sigma2.inv)
Sigma3 <- solve(Sigma3.inv)
Sigma4 <- solve(Sigma4.inv)
Sigma5 <- solve(Sigma5.inv)
Sigma6 <- solve(Sigma6.inv)
# # Sigma2 <- matrix(c(0.1,0,0,0.01),2,2)
# Sigma2 <- diag(c(0.1,0.01,0.1,0.01,0.1,0.01,0.1))
# # Sigma3 <- matrix(c(0.01,0,0,0.1),2,2)
# Sigma3 <- diag(c(0.01,0.1,0.01,0.1,0.01,0.1,0.01))

Example2 <- rbind(mvrnorm(n = 125, Mu1, Sigma1),
                  mvrnorm(n = 125, Mu2, Sigma2),
                  mvrnorm(n = 125, Mu3, Sigma3),
                  mvrnorm(n = 125, Mu4, Sigma4),
                  mvrnorm(n = 250, Mu5, Sigma5),
                  mvrnorm(n = 250, Mu6, Sigma6)
                  # data.unif,
                  # data.diamond
)

Example2 <- on.torus(Example2)

label <- c(rep(1,125),rep(2,125),
           rep(3,125),rep(4,125),
           rep(5,250),rep(6,250)
           # rep(4,50),
           # rep(5,50)
)

dat2 <- cbind( as.data.frame(Example2) , as.factor(label))
colnames(dat2) <- c("phi","psi","chi1","chi2","chi3", "chi4", "label")

Example2 <- rbind(mvrnorm(n = 125, Mu1, Sigma1),
                  mvrnorm(n = 125, Mu2, Sigma2),
                  mvrnorm(n = 125, Mu3, Sigma3),
                  mvrnorm(n = 125, Mu4, Sigma4),
                  mvrnorm(n = 250, Mu5, Sigma5),
                  mvrnorm(n = 250, Mu6, Sigma6)
                  # data.unif,
                  # data.diamond
)

Example2 <- on.torus(Example2)

label <- c(rep(1,125),rep(2,125),
           rep(3,125),rep(4,125),
           rep(5,250),rep(6,250)
           # rep(4,50),
           # rep(5,50)
)

dat2.test <- cbind( as.data.frame(Example2) , as.factor(label))
colnames(dat2.test) <- c("phi","psi","chi1","chi2","chi3","chi4","label")

pairs(dat2[, 1:6], pch = 19, col = rainbow(6)[label], upper.panel = NULL)


# testing codes---------------------------------------------------
type <- "he"

# predetermine J, alpha

data <- dat2[, -ncol(dat2)]
data.test <- dat2.test

data.test.label <- data.test[, ncol(dat2.test)]
data.test <- data.test[, -ncol(data.test)]

predicted.label <- matrix(NA, nrow = nrow(data), ncol = 2)

Jvec <- 10:50
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

cluster.num <- cluster.assign.number(data, Jmin = 10, Jmax = 50)
cluster.num$plot
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

Jhat <- 20
alphahat <- 0.1
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
# start <- Sys.time()
c <- cluster.assign.torus(data.test, icp.torus, level = alphahat, coord = t(combn(6,2)))
# end <- Sys.time()

# cat(end - start)

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

# GGally::ggpairs(result.dat[, 1:7], aes(color = result.dat[, 8]))
K <- length(unique(result.dat[,7]))
label <- as.factor(result.dat[, 7])
pairs(result.dat[, 1:4], pch = 19, col = rainbow(K)[label], upper.panel = NULL)

J <- 6

# case1 : extrinsic kmeans
kmeans.out <- KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
predicted.label1 <- predict_KMeans(cbind(cos(data.test),sin(data.test)) , kmeans.out$centroids)

adjustedRandIndex(predicted.label1,data.test.label)

# case2 : hierarchical
h.out <- hclust(d = ang.pdist(data), method = "complete")
predicted.label2 <- cutree(h.out, J)

adjustedRandIndex(predicted.label2,data.test.label)

# preparing for actual amino acid data ----------------------------
# high quality protein data obtained from PISCES server
# resolution : 1.6A(angstrom) or better
# R-factor : 0.22 or better
# Sequence percentage identity: <= 25%
protein <- read.table("PISCES_data.txt", header = TRUE)

# function for obtaining diheral angles of specified amino acid
# In Mardia (2012), the article uses ILE only data.
amino.data <- function(data, amino){
  IDs <- data$IDs
  ang <- data.frame()

  for (ID in IDs){
    id <- substr(ID, 1, 4)
    type <-substr(ID, 5, 5)

    tbl <- bio3d::torsion.pdb(bio3d::read.pdb(id))$tbl
    tbl <- data.frame(tbl, id = rownames(tbl)) %>%
      separate(id,into = c(',','position','type','rest'))

    tbl[,1:7] <- tbl[,1:7]/180*pi
    tbl[,1:7] <- on.torus(tbl[,1:7])
    tbl <- tbl %>% select(-",")

    tbl <- tbl[(tbl$type == type) & (tbl$rest == id), ]
    ang <- rbind(ang, tbl)
  }

  ang
}
load("ARG.Rdata")
load("ILE.Rdata")
sample.ILE <- sample(1:nrow(ILE), size = 2000)
ILE.sample <- ILE[sample.ILE, ]
pairs(ILE.sample)

sample.ARG <- sample(1:nrow(ARG), size = 1000)
ARG.sample <- ARG[sample.ARG, ]
pairs(ARG.sample)

type <- "ge"

Jvec <- 10:50
l <- list()

# sample spliting; preparing data
set.seed(2020)
n <- nrow(ILE.sample)
split.id <- rep(2,n)
split.id[ sample(n,floor(n/2)) ] <- 1

# testing for each J
for (j in Jvec){
  l[[j]] <- icp.torus.score(as.matrix(ILE.sample), split.id = split.id,
                            method = "kmeans",
                            kmeansfitmethod = type,
                            init = "k",
                            additional.condition = T,
                            param = list(J = j))
}

cluster.num <- cluster.assign.number(ILE.sample, Jmin = 10, Jmax = 50, method = type)
cluster.num$plot
cluster.num$cluster.number %>% group_by(Number) %>% count()

Jhat <- 25
alphahat <- 0.1
icp.torus <- l[[Jhat]]
# out <- rbind(out, data.frame(alpha = alphahat, J = Jhat))

c <- cluster.assign.torus(ILE.sample, icp.torus, level = alphahat, coord = c(1,2))

result.dat.out <- data.frame(ILE.sample, membership = as.factor(c$kmeans$cluster.id.outlier)) %>%
  mutate(membership = ifelse(membership == max(c$kmeans$cluster.id.outlier), "out", membership))

# GGally::ggpairs(result.dat[, 1:7], aes(color = result.dat[, 8]))
K <- length(unique(result.dat.out[, 5]))
label <- as.factor(result.dat.out[, 5])
pairs(result.dat.out[, 1:4], pch = 19, col = label, upper.panel = NULL)

result.dat.mahal <- data.frame(ILE.sample, membership = as.factor(c$kmeans$cluster.id.by.ehat)) %>%
  mutate(membership = ifelse(membership == max(c$kmeans$cluster.id.by.ehat), "out", membership))

K <- length(unique(result.dat.mahal[, 5]))
label <- as.factor(result.dat.mahal[, 5])
pairs(result.dat.mahal[, 1:4], pch = 19, col = label, upper.panel = NULL)
