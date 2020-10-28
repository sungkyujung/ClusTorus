
torsion_to_dat <- function(pdb.6M16){
tbl <- data.frame(phi = pdb.6M16$phi, psi = pdb.6M16$psi, id = rownames(pdb.6M16$tbl)) %>% 
  separate(id,into = c(',','position','type','rest'))
tbl
}

plot_clustering <- function(icp.torus, choosealpha){
  
  ia <- icp.torus.eval(icp.torus, level = choosealpha, eval.point = grid.torus())
  c <- cluster.assign.torus(data, icp.torus, level = choosealpha) 
  
  b <- data.frame(ia$phi,ia$psi, ia$Chat_mix == 1, ia$Chat_max == 1, ia$Chat_e == 1)
  colnames(b) <- c("phi","psi","C_mix","C_max","C_e")
  b<- b %>%  pivot_longer(3:5, names_to = "Type", values_to = "value") %>% filter(Type == "C_e") 
  g0 <- ggplot() + 
    geom_contour(aes(phi, psi, z = ifelse(value,1,0)), data = b, size = 1,lineend = "round" ) + 
    geom_point(mapping = aes(x,y,shape = c, color = c), data = data.frame(x = data[,1],y =data[,2], 
                                                                          c = as.factor(c$cluster.id.by.ehat))) 
  g0 + ggtitle(paste("ICP, vM2 mixture with J =",chooseJ,", alpha =",choosealpha)) 
  
  g2 <- g0
  
  level = choosealpha
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
  g2 <- g2 + 
    scale_x_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","2pi/3","2pi"), limits = c(-eps,2*pi+eps))+ 
    scale_y_continuous(breaks = c(0,1,2,3,4)*pi/2, labels = c("0","pi/2","pi","2pi/3","2pi"), limits = c(-eps,2*pi+eps))+
    ggtitle(paste("alpha=",choosealpha)) 
  
  g2
}
