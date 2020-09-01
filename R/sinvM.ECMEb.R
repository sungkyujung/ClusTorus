#'
sinvM.ECMEb <- function(wstat, kappa1, kappa2, lambda, mu1, mu2, THRESHOLD = 1e-10){

  mu1 <- c(cos(mu1),sin(mu1))
  mu2 <- c(cos(mu2),sin(mu2))

  r <- 0
  while(TRUE){
    mu1c <- kappa1 * wstat$y1bar + lambda * wstat$S12 %*% mu2
    mu1c <- mu1c / sqrt(sum(mu1c^2))
    mu2c <- kappa2 * wstat$y2bar + lambda * t(wstat$S12) %*% mu1
    mu2c <- mu2c / sqrt(sum(mu2c^2))

    # STOP IF NOT CHANGING
    diffb <- norm(mu1c - mu1) + norm(mu2c - mu2)
    r <- r + 1
    if(r >= 100 | diffb < THRESHOLD){break}
    mu1 <- mu1c
    mu2 <- mu2c

  }
  mu1 <- atan2(mu1[2], mu1[1])
  mu2 <- atan2(mu2[2], mu2[1])
  return(list(mu1 = ifelse(mu1 < 0, mu1 + 2*pi, mu1),
              mu2 = ifelse(mu2 < 0, mu2 + 2*pi, mu2)))
}
