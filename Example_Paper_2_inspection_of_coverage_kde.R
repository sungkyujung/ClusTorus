library(MASS)
library(cowplot)
library(tidyverse)
library(ClusTorus)
set.seed(2)
Out.dt <- NULL
RR  = 100
for (i in 1:RR){
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


# Evaluate Coverage -------------------------------------------------------


  data <- dat1[,1:2]
  data.test <- dat1.test
  data.test.label <- data.test[,3]
  data.test <- as.matrix(data.test[,1:2])

  grid.test <- grid.torus()
  testing.n <- nrow(data)


  n2 <- floor(nrow(data)/2)
  alphavec <- 1:floor(n2/2) / n2
  alphavec <- alphavec[alphavec <= 0.4]
  N <- length(alphavec)



  icp.torus<- icp.torus.score(as.matrix(data), split.id = NULL,
                                  method = "a",
                                  mixturefitmethod = "a",
                                  kmeansfitmethod = "g",
                                  param = list(J = 13, concentration = 25))
  CC <- icp.torus.eval(icp.torus, level = alphavec, eval.point = grid.torus())



  C <- CC$Chat_e
  Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  for (j in 1:testing.n){
   Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  }
  Out.dt <- rbind(Out.dt, data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_e"))

  C <- CC$Chat_kde
  Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  for (j in 1:testing.n){
   Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  }
  Out.dt <- rbind(Out.dt,  data.frame(alpha = alphavec,coverage = colMeans(Inclusion.test), Type = "C_KDE")
  )


  C <- CC$Chat_mix
  Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  for (j in 1:testing.n){
   Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  }
  Out.dt <- rbind(Out.dt,
                 data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_mix"))

  C <- CC$Chat_max
  Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  for (j in 1:testing.n){
    Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  }
  Out.dt <- rbind(Out.dt,
                  data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_max"))

  C <- CC$Chat_kmeans
  Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
  for (j in 1:testing.n){
    Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
  }
  Out.dt <- rbind(Out.dt,
                  data.frame(alpha = alphavec,coverage = colMeans(Inclusion.test), Type = "C_kmeans"))
}

Out.dt %>% arrange(alpha,Type) %>%  head()

CIdat <- data.frame(alpha = alphavec, coverage = (1 - alphavec) - 2*sqrt(alphavec*(1-alphavec)/testing.n), Type = "CI(95%)")

g1 <- Out.dt %>%
  group_by(alpha,Type) %>%
  summarise(n = n(), coverage_mean = mean(coverage), coverage_sd = sd(coverage)) %>%
  ggplot(aes(x = 1-alpha, y = coverage_mean, color = Type)) +
  geom_line(aes(1-alpha, coverage), data = CIdat, color = "gray", size = 2) +
  geom_errorbar(aes(ymin = coverage_mean - 2*coverage_sd, ymax = coverage_mean + 2*coverage_sd) )+
  geom_abline(slope = 1, intercept = 0) +
  geom_line() + ggtitle("Model I")

g1

Out.dt1 <- Out.dt
################



Out.dt <- NULL

for (i in 1:RR){
# Prepare data 2 ------------------------------------------------------------

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


data <- dat2[,1:2]
data.test <- dat2.test

data.test.label <- data.test[,3]
data.test <- as.matrix(data.test[,1:2])


grid.test <- grid.torus()
testing.n <- nrow(data)

n2 <- floor(nrow(data)/2)
alphavec <- 1:floor(n2/2) / n2
alphavec <- alphavec[alphavec <= 0.4]
N <- length(alphavec)

icp.torus<- icp.torus.score(as.matrix(data), split.id = NULL,
                            method = "all",
                            mixturefitmethod = "a",
                            kmeansfitmethod = "g",
                            param = list(concentration = 25, J = 13))
CC <- icp.torus.eval(icp.torus, level = alphavec, eval.point = grid.torus())


C <- CC$Chat_e
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt, data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_e"))

C <- CC$Chat_kde
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt,  data.frame(alpha = alphavec,coverage = colMeans(Inclusion.test), Type = "C_KDE"))


C <- CC$Chat_mix
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt,
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_mix"))

C <- CC$Chat_max
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt,
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_max"))


C <- CC$Chat_kmeans
Inclusion.test <- matrix(0, nrow = testing.n, ncol = ncol(C))
for (j in 1:testing.n){
  Inclusion.test[j,] <- C[which.min(rowSums( (t( t(grid.test) - data.test[j,]))^2) ),]
}
Out.dt <- rbind(Out.dt,
                data.frame(alpha = alphavec, coverage = colMeans(Inclusion.test), Type = "C_kmeans"))
}

Out.dt2 <- Out.dt

Out.dt %>% arrange(alpha,Type) %>%  head()

CIdat <- data.frame(alpha = alphavec, coverage = (1 - alphavec) - 2*sqrt(alphavec*(1-alphavec)/testing.n), Type = "CI(95%)")

g2 <- Out.dt %>%
  group_by(alpha,Type) %>%
  summarise(n = n(), coverage_mean = mean(coverage), coverage_sd = sd(coverage)) %>%
  ggplot(aes(x = 1-alpha, y = coverage_mean, color = Type)) +
  geom_line(aes(1-alpha, coverage), data = CIdat, color = "gray", size = 2) +
  geom_errorbar(aes(ymin = coverage_mean - 2*coverage_sd, ymax = coverage_mean + 2*coverage_sd) )+
  geom_abline(slope = 1, intercept = 0) +
  geom_line() + ggtitle("Model II")
g2

#
# g2 <- Out.dt %>%
#   group_by(alpha,Type) %>%
#   summarise(n = n(), coverage_mean = mean(coverage), coverage_sd = sd(coverage)) %>%
#   ggplot(aes(x = 1-alpha, y = coverage_mean, color = Type)) +
#   geom_errorbar(aes(ymin = coverage_mean - 2*coverage_sd, ymax = coverage_mean + 2*coverage_sd) )+
#   geom_abline(slope = 1, intercept = 0) +
#   geom_line() + ggtitle("Model II")

save.image(file = "Example_Paper_2_simulation_data2.RData")

load(file = "Example_Paper_2_simulation_data2.RData")
plot_grid(g1 ,
          g2 , label_size = 12)
ggsave("./examples/Toy_Data4.png", width = 8, height = 4)


# variance decomposition --------------------------------------------------
n1 <- 270

Out.dt1 %>%
  group_by(alpha,Type) %>%
  summarise(n = n(),
            coverage_mean = mean(coverage),
            coverage_var = var(coverage)) %>%
  mutate(var_test = alpha *(1-alpha)/n1,
         var_train = coverage_var - var_test) %>%
  mutate(Phat = var_test / coverage_var*100,
         Cn = var_train / coverage_var * 100) %>%
  select(-coverage_mean,-var_train,-var_test) %>%
  pivot_longer(cols = 5:6, names_to = "Source", values_to = "Variance") %>%
  ggplot(aes(x= 1-alpha, y = Variance, fill = Source)) + geom_area() + facet_grid(Type ~.)


# variance decomposition --------------------------------------------------
n1 <- 500

Out.dt2 %>%
  group_by(alpha,Type) %>%
  summarise(n = n(),
            coverage_mean = mean(coverage),
            coverage_var = var(coverage)) %>%
  mutate(var_test = alpha *(1-alpha)/n1,
         var_train = coverage_var - var_test) %>%
  mutate(Phat = var_test / coverage_var*100,
         Cn = var_train / coverage_var * 100) %>%
  dplyr::select(-coverage_mean,-var_train,-var_test) %>%
  pivot_longer(cols = 5:6, names_to = "Source", values_to = "Variance") %>%
  ggplot(aes(x= 1-alpha, y = Variance, fill = Source)) + geom_area() + facet_grid(Type ~.)


Out_all <- rbind( cbind(Out.dt1, Model = "I", n1 = 270),
                  cbind(Out.dt2, Model = "II", n1 = 500))

Out_all %>%
  group_by(Model,alpha,Type) %>%
  summarise(n = n(),
            coverage_mean = mean(coverage),
            coverage_var = var(coverage), n1 = min(n1)) %>%
  mutate(var_test = alpha *(1-alpha)/n1,
         var_train = coverage_var - var_test) %>%
  mutate(Phat = var_test / coverage_var*100,
         Cn = var_train / coverage_var * 100) %>%
  dplyr::select(-coverage_mean,-var_train,-var_test,-n1) %>%
  pivot_longer(cols = 6:7, names_to = "Source", values_to = "Variance") %>%
  ggplot(aes(x= 1-alpha, y = Variance, fill = Source)) + geom_area() + facet_grid(Type ~ Model)

ggsave("./examples/Toy_Data5.png", width = 6, height = 8)

Out_all %>%
  group_by(Model,alpha,Type) %>%
  summarise(n = n(),
            coverage_mean = mean(coverage),
            coverage_var = var(coverage), n1 = min(n1)) %>%
  mutate(var_test = alpha *(1-alpha)/n1,
         var_train = coverage_var - var_test) %>%
  mutate(VPhat = var_test,
         VCn = var_train) %>%
  dplyr::select(-coverage_mean,-var_train,-var_test,-n1) %>%
  pivot_longer(cols = 6:7, names_to = "Source", values_to = "Variance") %>%
  ggplot(aes(x= 1-alpha, y = Variance, fill = Source)) + geom_area() + facet_grid(Type ~ Model)

ggsave("./examples/Toy_Data6.png", width = 6, height = 8)

