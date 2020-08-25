library(tidyverse)
library(bio3d)
source('routines.R')
source('routines_MISC.R')
library(cowplot)
# bio data ------------------------------------------ 
# pdb.6M15 <- torsion.pdb(read.pdb("6M15"))  # HKU2
# pdb.6M16 <- torsion.pdb(read.pdb("6M16") ) # SADS-Cov
# pdb.6VXX <- torsion.pdb(read.pdb("6VXX"))  # SARS-Cov2 (Covid-19)
# save(pdb.6M15,pdb.6M16,pdb.6VXX,file = "coronavirusdata.RData") 
load(file = "coronavirusdata.RData")
head(pdb.6VXX$tbl)
dataXX <- torsion_to_dat(pdb.6VXX)
dataXX %>% group_by(type) %>% summarize(n = n(),maxposition = max(position))
dataXX %>% filter(type %in% c("A","B")) %>% group_by(type) %>% summarize(n = n(),maxposition = max(position))
dataXX <- dataXX %>% filter(type %in% c("B"))
head(dataXX) 
tail(dataXX)
plot(dataXX$position)

n <- nrow(dataXX)
n1 <- n/2
n2 <- n/2

dataXX %>% ggplot(aes(phi, psi)) + geom_point(size = 1) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*90, labels = c("-pi","-pi/2","0","pi/2","pi"))+ 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*90, labels = c("-pi","-pi/2","0","pi/2","pi"))  
ggsave("./coronavirus/RamPlotXX.png", width = 6, height = 4)


data16 <- torsion_to_dat(pdb.6M16)
head(data16) 
data16 <- data16 %>% filter(type == "B") %>% select(phi,psi, position)

data15 <- torsion_to_dat(pdb.6M15)
head(data15) 
data15 <- data15 %>% filter(type == "B") %>% select(phi,psi, position)
nrow(data15)

dataXX %>% ggplot(aes(phi, psi)) + geom_point(size = 1) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2)*90, labels = c("-pi","-pi/2","0","pi/2","pi"))+ 
  scale_y_continuous(breaks = c(-2,-1,0,1,2)*90, labels = c("-pi","-pi/2","0","pi/2","pi")) + 
  geom_point(aes(phi,psi), data = data15, color = "blue",size = 0.5)+ 
  geom_point(aes(phi,psi), data = data16, color = "red",size = 0.5)
ggsave("./coronavirus/RamPlot_XX.png", width = 6, height = 4)




data <- data.frame(phi = dataXX$phi/180*pi, psi = dataXX$psi/180*pi)
# data <- data[-which(is.na(data[,1])|is.na(data[,2])),]
split.id <- rep(c(1,2),n/2)
data <- on.torus(as.matrix(data[,1:2]))

plot(data)
nrow(data)

# Clustering by existing methods ------------------------------------------
J <- 3
filename <- "./coronavirus/6VXX_"
# source("Ex_Master.R") -------------------

# Perfomr naive clustering 
library(tidyverse)
library(ClusterR)
library(mclust)  
## By naive K-means
kmeans.out<-KMeans_rcpp(data,clusters = J)
g_kmeans <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(kmeans.out$clusters)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle("K-means (ignoring the angular constraint)")
g_kmeans
ggsave(paste(filename,"1kmeans.png",sep = ""), width = 6, height = 4)


## By Non-angular Gaussian Mixture (using state-of-the-art mclust)
# BIC <- mclustBIC(data, G = J)
# plot(BIC)
# mod1 <- Mclust(data, x = BIC)
# summary(mod1, parameters = TRUE)
# plot(mod1, what = "classification",xlab = "phi",ylab = 'psi', main = "Normal mixture") 

# define angular distance: 

pdist.data2 <- ang.pdist(data) # Use the pairwise L2-angular distance for clustering 

## PAM (Partitioning around medoids - Kaufman and Rousseeuw(1990) )
pam.out <- Cluster_Medoids(as.matrix(pdist.data2),clusters = J)
g_pam <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(pam.out$clusters)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle("Partitioning around medoids")
g_pam 
ggsave(paste(filename,"2pam.png",sep = ""), width = 6, height = 4)

## K-means in the ambient space with kmeans++ initialization

kmeans.out<-KMeans_rcpp(cbind(cos(data),sin(data)),clusters = J)
g_kmeans2 <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(kmeans.out$clusters)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle("K-means with chordal distance (in the ambient space) ")
g_kmeans2
ggsave(paste(filename,"2kmeans_amb.png",sep = ""), width = 6, height = 4)

## Hierarchical
hc.complete <- hclust(pdist.data2, method="complete")
library("ggdendro")
# ggdendrogram(hc.complete, rotate=TRUE, size=2) + labs(title="Complete Linkage")
#ggdendrogram(hc.complete,main="Average Linkage", xlab="", sub="", cex=.9)
membership <- cutree(hc.complete, J)
g_hier <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(membership)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle("Hierachical clustering with average L2-Angular distance")
g_hier
ggsave(paste(filename,"3hier.png",sep = ""), width = 6, height = 4)

eps <- 0.0

plot_grid(g_kmeans + 
            scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
            scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))
          , g_pam + 
            scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
            scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))
          , g_kmeans2 + ggtitle("K-means (in the ambient space)") + 
            scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
            scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))
          , g_hier  + ggtitle("Hierachical clustering")  + 
            scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
            scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))
          , label_size = 12)

ggsave(paste(filename,"11figure_clust.png",sep = ""), width = 8, height = 4*8/6)


# Find alpha and J --------------------------------------------------------

Jvec <- 3:65
l <- list()
set.seed(0)
n <- nrow(data)
n2 <- sum(split.id == 2)
 for (j in Jvec){
   l[[j]] <- icp.torus.score(data, split.id = split.id,
                             method = "mixture",
                             mixturefitmethod = "a",
                             param = list(J = j))
 }
alphavec <- 1:floor(n2/2) / n2
N <- length(alphavec)

# need a data frame (alpha, J, mu, alpha + mu)
 out <- data.frame()
 for (j in Jvec){
   Mvec <- alphavec
   a<-icp.torus.eval(l[[j]], level = alphavec, eval.point = grid.torus())
   for (i in 1:N){
     Mvec[i] <- sum(a$Chat_e[,i])/10000
   }
   out <- rbind(out, data.frame(alpha = alphavec, J = j, mu = Mvec, criterion = alphavec + Mvec))
 } 
 save(l,out,file = "coronavirus_mix_fit.RData")

load(file = "coronavirus_mix_fit.RData")
out %>% ggplot(aes(x= alpha, y = mu, color = J)) + geom_point()
ggsave(paste(filename,"7icp_areaC.png",sep = ""), width = 6, height = 4)

out.index <- which.min(out$criterion)
out[out.index,] 

Jhat <- out[out.index,2]
alphahat <- out[out.index,1]

out %>% filter(alpha == alphahat) %>% ggplot(aes(J,mu)) + geom_point() + geom_line()

# ICP-Clustering by chosen alpha and J ------------------------------------
icp.torus <- l[[Jhat]]

temp <- data.frame(iteration = 1:length(icp.torus$mixture$fit$loglkhd.seq), LKHD = icp.torus$mixture$fit$loglkhd.seq, Delta = c(NA, (  rowSums(diff(icp.torus$mixture$fit$param.seq)^2)))) 

temp %>% ggplot(aes(iteration, LKHD)) + geom_point() + geom_line() + 
  ggtitle(paste("vM2 Mixture fit with J =",Jhat,"components")) 
ggsave(paste(filename,"8mixture_fit1.png",sep = ""), width = 6, height = 4)

temp %>% ggplot(aes(iteration,Delta)) +
  geom_point() + geom_line() +
  scale_y_log10() + 
  ggtitle(paste("vM2 Mixture fit with J =",Jhat,"components")) 
ggsave(paste(filename,"8mixture_fit2.png",sep = ""), width = 6, height = 4)


ia <- icp.torus.eval(icp.torus, level = alphahat, eval.point = grid.torus())
b <- data.frame(ia$phi,ia$psi, ia$Chat_mix == 1, ia$Chat_max == 1, ia$Chat_e == 1)
colnames(b) <- c("phi","psi","C_mix","C_max","C_e")
head(b)

b<- b %>%  pivot_longer(3:5, names_to = "Type", values_to = "value")
g0 <- ggplot() + geom_contour(aes(phi, psi, z = ifelse(value,1,0),linetype = Type, color =Type), data = b, size = 1,lineend = "round" ) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2])) + 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
g0 + ggtitle(paste("ICP, vM2 mixture with J=",Jhat,", alpha=",alphahat)) 
ggsave(paste(filename,"9icp.png",sep = ""), width = 6, height = 4)


c <- cluster.assign.torus(data, icp.torus, level = alphahat) 
g_e <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.ehat)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle(paste("Clustering by hat{e}, K=",c$ncluster)) + 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
g_e
ggsave(paste(filename,"10clustering1.png",sep = ""), width = 6, height = 4)

membership <- c$cluster.id.outlier 
membership <- ifelse(membership == max(unique(c$cluster.id.outlier)),"out",membership)

g1 <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(membership)) %>%  
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle(paste("Clusters and outliers, K=",length(unique(c$cluster.id.outlier)))) + 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))
g1

ggsave(paste(filename,"10clustering2.png",sep = ""), width = 6, height = 4)


g_mah <- data.frame(phi = data[,1], psi = data[,2], membership = as.factor(c$cluster.id.by.Mah.dist)) %>% 
  ggplot(aes(phi,psi, color = membership)) + geom_point() + ggtitle(paste("Clustering by Mahalanobis distance, K=",c$ncluster))
g_mah  

ggsave(paste(filename,"10clustering3.png",sep = ""), width = 6, height = 4)


g2 <- g_e

level = alphahat
n2 <- icp.torus$n2
ialpha <- floor( (n2 + 1) * level)
t <- icp.torus$mixture$score_ellipse[ialpha]
# Draw.ellipses.bvn.approx(data, icp.torus$mixture$fit$parammat, t, data, c$cluster.id.outlier)

ellipse.param <- icp.torus$mixture$ellipsefit
J <- length(ellipse.param$mu1)

# all_nine_ellipses
theta <- seq(0,2*pi,length.out = 999)
Z <- cbind(cos(theta), sin(theta))

shift <- matrix(0,ncol = 2, nrow = 9)
shift[,1] <- c(0,2*pi,-2*pi)
shift[,2] <- rep(c(0,2*pi,-2*pi), each = 3)


for(j in 1:J){
  mu <- c(ellipse.param$mu1[j], ellipse.param$mu2[j])
  Sinv <- ellipse.param$Sigmainv[[j]]
  c.minus.t <- ellipse.param$c[j] - t
  
  if(c.minus.t < 0){
    cat("skip",j,",")
    next}
  cat("draw",j,",")
  M <- eigen(Sinv/c.minus.t)
  Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
  R <- Mmhalf %*% t(Z) 
  for( shift.id in 1:9){
    RR <- R + mu + shift[shift.id,]
    g2 <-   g2 + geom_polygon(aes(x = phi, y = psi),color = "blue",alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
  }
  
}

eps <- pi/10
g2 + 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))

ggsave(paste(filename,"10clustering4.png",sep = ""), width = 6, height = 4)



# a figure for paper ------------------------------------------------------


plot_grid(g_cp, g0, g_e + 
            scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
            scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))
          , g1 + 
            scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))+ 
            scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(-eps,2*pi+eps))
          , label_size = 12)

ggsave(paste(filename,"11figure_showcase1.png",sep = ""), width = 6, height = 4)
ggsave(paste(filename,"11figure_showcase2.png",sep = ""), width = 8, height = 4*8/6)


####################


out[out.index,]
Jhat
alphahat
library(cowplot)


#icp.torus_general <- icp.torus.score(data, split.id = NULL,
#                            method = "mixture",
#                            mixturefitmethod = "g",
#                            param = list(J = 10))

#source('plot_clustering.R')
# p1 <- plot_clustering(icp.torus_general, choosealpha = 0.02)
# p2 <- plot_clustering(icp.torus_general, choosealpha = 0.05)
# p3 <- plot_clustering(icp.torus_general, choosealpha = 0.10)
# p4 <- plot_clustering(icp.torus_general, choosealpha = 0.20)
# 
# plot_grid(p1, p2, p3,p4, label_size = 12)
# ggsave(paste(filename,"6M16_ClusterPlot_general.png",sep=""), width = 12, height = 8)
# 


chooseJ <- Jhat 
icp.torus <- l[[chooseJ]]
source('plot_clustering.R')
p1 <- plot_clustering(icp.torus, choosealpha = 0.02)
p2 <- plot_clustering(icp.torus, choosealpha = 0.03)
p3 <- plot_clustering(icp.torus, choosealpha = 0.05)
p4 <- plot_clustering(icp.torus, choosealpha = 0.10)

plot_grid(p1, p2, p3,p4, label_size = 12)
ggsave(paste(filename,"ClusterPlot.png",sep=""), width = 12, height = 8)

# kappa and alpha 
# Now with kde-based fit

set.seed(101)
n <- nrow(data)

icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 25))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.02, eval.point = grid.torus())
p1<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2]))+ 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  ggtitle('KDE: k = 25, alpha = 0.02')

icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.1, eval.point = grid.torus())
p2<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2]))+ 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  ggtitle('KDE: k = 25, alpha = 0.1')


icp.torus.kde<- icp.torus.score(data, split.id = split.id,
                                method = "kde",
                                mixturefitmethod = "a",
                                param = list(concentration = 100))
icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.02, eval.point = grid.torus())
p3<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2]))+ 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  ggtitle('KDE: k = 100, alpha = 0.02')

icp.kde.region <- icp.torus.eval(icp.torus.kde, level = 0.1, eval.point = grid.torus())
p4<- data.frame(Chat = as.vector(icp.kde.region$Chat_kde), psi = icp.kde.region$psi, phi = icp.kde.region$phi) %>% 
  ggplot() + geom_contour(aes(phi,psi,z = Chat)) + 
  geom_point(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2]))+ 
  scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+ 
  scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","3pi/2","2pi"), limits = c(0,2*pi))+
  ggtitle('KDE: k = 100, alpha = 0.1')

plot_grid(p1, p2, p3,p4, label_size = 12)
ggsave(paste(filename,"ICP_KDE.png",sep=""), width = 12, height = 8)



# choice of alpha and j browse options ------------------------------------
concentration.vec <- 10^( seq(log10(4),log10(150),length.out = 25) )
nkappa <- length(concentration.vec)
mu.vec1 <- concentration.vec
mu.vec2 <- concentration.vec

n <- nrow(data)

for (j in 1:nkappa){
  icp.torus.kde<- icp.torus.score(data,method = "kde",split.id = split.id ,param = list(concentration = concentration.vec[j]))
  icp.kde.region <- icp.torus.eval(icp.torus.kde, level = c(0.05,0.1), eval.point = grid.torus()) 
  mu.vec1[j] <- sum(icp.kde.region$Chat_kde[,1])/10000
  mu.vec2[j] <- sum(icp.kde.region$Chat_kde[,2])/10000
}

pICP_KDE_k <- data.frame(mu1 = mu.vec1, mu2 = mu.vec2, concentration = concentration.vec) %>% 
  pivot_longer(1:2, names_to = "Level") %>% mutate(Level = as.factor ( ifelse(Level=="mu1",0.05,0.1)) ,mu = value) %>% 
  ggplot(aes(concentration, mu, color = Level, group = Level)) + 
  geom_point() + geom_line()  + scale_x_log10()
pICP_KDE_k
ggsave(paste(filename,"ICP_tuning_1.png",sep=""), width = 6, height = 4)


# alpha vs j for MIXTURE --------------------------------------------------


out[out.index,]

pICP_MIX_J <- out %>% filter( abs(alpha - 0.1) < 0.001 | abs(alpha - 0.05) < 0.001 ) %>% 
  mutate(Level = as.factor(ifelse( abs(alpha - 0.1) < 0.001, 0.1, 
                                   0.05))) %>%  ggplot(aes(J,mu, color = Level)) + geom_point() + geom_line()
pICP_MIX_J
ggsave(paste(filename,"ICP_tuning_2.png",sep=""), width = 6, height = 4)

plot_grid(pICP_KDE_k + ggtitle('KDE, over kappa'), 
          pICP_MIX_J+ ggtitle('Mixture, over J'), label_size = 12)
ggsave(paste(filename,"ICP_tuning_3.png",sep=""), width = 6, height = 3)



out %>% group_by(alpha) %>% ggplot(aes(alpha,criterion,color = J)) + 
  geom_point() + scale_colour_distiller(palette="Spectral") +geom_abline(slope = 0, intercept = out[out.index,4])

ggsave(paste(filename,"ICP_tuning_4criterion.png",sep=""), width = 6, height = 3)

out %>%  ggplot(aes(alpha,1-mu,color = J, group = J)) + geom_line() +
  geom_point() + scale_colour_distiller(palette="Spectral") +geom_abline(slope = 1, intercept = 1-out[out.index,4])
ggsave(paste(filename,"ICP_tuning_4criterion_ud.png",sep=""), width = 6, height = 3)


out %>% group_by(J) %>% summarise(AUC = mean(1-mu)/2) %>% 
  ggplot(aes(J,AUC,color = J)) + geom_line()+geom_point() + 
  scale_color_gradientn(colours = rainbow(length(Jvec)))

out %>% group_by(J) %>% summarise(AUC = mean(1-mu)/2) %>% arrange(desc(AUC))


out %>%  ggplot(aes(alpha,1-mu,color = J, group = J)) + geom_line() +
  geom_point() + scale_colour_distiller(palette="Spectral") +geom_abline(slope = 1, intercept = 1-out[out.index,4]) +
  geom_line(aes(alpha,1-mu),data = out %>% filter(J == 19), color = "black", size = 2)


out %>%  ggplot(aes(alpha,1-mu,color = J, group = J)) + geom_line() +
  geom_point() + scale_color_gradientn(colours = rainbow(length(Jvec)))  +geom_abline(slope = 1, intercept = 1-out[out.index,4]) +
  geom_line(aes(alpha,1-mu),data = out %>% filter(J == 19), color = "black", size = 2)


out %>% filter(J == 19) %>% ggplot(aes(alpha,mu)) + 
  geom_point() + geom_abline(slope = -1, intercept = out[out.index,4])


# inspect a little further on out 

best.option <- out %>% as.data.frame()  %>% arrange(criterion) %>% head(1) 
best2.option <- out %>% as.data.frame()  %>% filter(J != Jhat) %>% arrange(criterion) %>% head(1) 
best.option
best2.option

icp.torus <- l[[best.option$J]]
plot_clustering(icp.torus, choosealpha = best.option$alpha)



icp.torus <- l[[best2.option$J]]
plot_clustering(icp.torus, choosealpha = best2.option$alpha)

ggsave(paste(filename,"12figure_alternate_J18alpha087.png",sep = ""), width = 8, height = 4*8/6)







# 
# ggplot() + geom_point(mapping = aes(x,y,shape = c, color = c), data = data.frame(x = data[,1],y =data[,2], 
#                                                                                  c = as.factor(c$cluster.id.by.ehat))) +
#   geom_path(mapping = aes(x,y), data = data.frame(x = data[,1],y =data[,2], 
#                                                   c = as.factor(c$cluster.id.by.ehat))) 
# 
# data.frame(id = 1:nrow(data), x = data[,1],y =data[,2], 
#            c = (c$cluster.id.by.ehat)) %>% ggplot(aes(id,c,color = c)) + geom_point()  + geom_line()+
# geom_line(aes(id,c), 
#           data = 
#             data.frame(id = 1:nrow(data), x = data[,1],y =data[,2], 
#                        c = (c$cluster.id.by.ehat)))

cc <- c$cluster.id.by.ehat 
cMat <- matrix(0, nrow= c$ncluster, ncol = c$ncluster)
for( r in 1:(length(cc)-1)){
  cMat[cc[r], cc[r+1]] = cMat[cc[r], cc[r+1]] + 1 
}
cMat
cMat / rowSums(cMat)
length(cc)
 

# rownames(cMat) <- c("alpha-helix(r)","beta-sheet","alpha-helix(l)","Misc.")
# colnames(cMat) <- c("alpha-helix(r)","beta-sheet","alpha-helix(l)","Misc.")
# library("qgraph")
# Graph_cor <- qgraph(cMat,label.prop=1.3,layout="spring",theme='TeamFortress',aspect=0.3)


#graph_from_adjacency_matrix()
#cMatg<-graph_from_adjacency_matrix(cMat, mode = "directed", weighted = TRUE,
#                            diag = TRUE, add.colnames = NULL, add.rownames = NA)
#E(cMatg)$weight
#plot.igraph(cMatg, edge.arrow.size = 0.2)
#plot(cMatg)
 