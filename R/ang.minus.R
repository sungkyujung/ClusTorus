#' Augular subtraction
#'
#' \code{ang.minus()} computes angular subtraction.
#'
#' @param x,y angular data on \eqn{[0, 2\pi)}
#' @export
#' @examples
#' \dontrun{
#' x <- pi/2
#' y <- pi/3
#'
#' ang.minus(x, y)
#' }
ang.minus <- function(x,y){
  # x "minus" y if x and y are [0, 2*pi)
  t <- rbind(x - y, x - y + 2 * pi, x - y - 2 * pi)

  tind <- apply(abs(t), 2, which.min)

  tt <- t[1, ]
  tt[tind == 2] <- t[2, tind == 2]
  tt[tind == 3] <- t[3, tind == 3]

  tt
}
