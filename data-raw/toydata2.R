## prepare `toydata2` dataset goes here
library(MASS)
library(ClusTorus)

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
toydata2 <- rbind(dat2, dat2.test)
colnames(toydata2) <- c("phi","psi","label")

usethis::use_data(toydata2, overwrite = TRUE)
