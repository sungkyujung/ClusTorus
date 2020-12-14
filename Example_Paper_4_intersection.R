library(latex2exp)

g2 <- ggplot()

mu1 = c(0.1,0.1)*2*pi
mu2 = c(0.8,0.2)*2*pi
Sigmainv1 <- matrix(c(1,0,
                      0,1),
                        nrow = 2,ncol = 2)
Sigmainv2 <- 0.5*matrix(c(2,1,
                      1,2),
                    nrow = 2,ncol = 2)

theta <- seq(0,2*pi,length.out = 999)
Z <- cbind(cos(theta), sin(theta))
shift <- matrix(0,ncol = 2, nrow = 9)
shift[,1] <- c(0,2*pi,-2*pi)
shift[,2] <- rep(c(0,2*pi,-2*pi), each = 3)

# E1:
mu <- mu1
Sinv <- Sigmainv1
c.minus.t <- 1
M <- eigen(Sinv/c.minus.t)
Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
R <- Mmhalf %*% t(Z)
for( shift.id in 1:9){
  RR <- R + mu + shift[shift.id,]
  g2 <-   g2 + geom_polygon(aes(x = phi, y = psi),color = "blue",alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
}

# E2:
mu <- mu2
Sinv <- Sigmainv2
c.minus.t <- 1
M <- eigen(Sinv/c.minus.t)
Mmhalf <- M$vectors %*% diag( sqrt(1/M$values) ) %*% t(M$vectors)
R <- Mmhalf %*% t(Z)
RR <- R + mu
g2 <-   g2 + geom_polygon(aes(x = phi, y = psi),color = "red",alpha = 0.1, data = data.frame(phi = RR[1,],psi = RR[2,], value = 1))
g2

g2 + geom_vline(xintercept = c(-1,0,1,2)*2*pi) +
     geom_hline(yintercept = c(-1,0,1,2)*2*pi) +
  scale_x_continuous(breaks = c(-1,0,1,2)*2*pi,
                     labels = c("-2pi","0","2pi",""),
                     limits = c(-2.2*pi,4*pi)) +
  scale_y_continuous(breaks = c(-1,0,1,2)*2*pi,
                     labels = c("-2pi","0","2pi",""),
                     limits = c(-2.2*pi,4*pi)) +
  annotate(geom='text', x=mu1[1], y=mu1[2],
           label=latex2exp("$\\tilde{E}_1", output='character'), parse=TRUE) +
  annotate(geom='text', x=mu2[1], y=mu2[2],
           label=latex2exp("$\\tilde{E}_2", output='character'), parse=TRUE) +
  annotate(geom='text', x=mu1[1]+2*pi, y=mu1[2],
           label=latex2exp("$\\tilde{E}_1'", output='character'), parse=TRUE)

ggsave("./examples/Toy_Figure_intersection2.png", width = 6, height = 6)
ggsave("./examples/Toy_Figure_intersection.png", width = 4, height = 4)
