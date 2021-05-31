## prepare `toydata1` dataset goes here
library(MASS)
library(ClusTorus)
set.seed(20201)

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
toydata1 <- rbind(dat1, dat1.test)
colnames(toydata1) <- c("phi","psi","label")

usethis::use_data(toydata1, overwrite = TRUE)
