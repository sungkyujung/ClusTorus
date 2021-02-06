#' Grid on torus
#'
#' \code{grid.torus} returns an equally-spaced grid on torus.
#'
#' @param d number for dimension. Default is 2.
#' @param grid.size number of grid for each axis. Default value is 100.
#' @return returns (grid.size) x (grid.size) numeric matrix
#'   which indicates the grid points on torus.
#' @export
#' @examples
#' grid.torus(d = 2, grid.size = 100)
grid.torus <- function(d = 2 , grid.size = 100){
  # returns grid points on torus of size (grid.size) x (grid.size)

  grid <- matrix(0, ncol = d, nrow = grid.size^d)

  Axis <- seq(0, 2 * pi, length = grid.size)
  for (i in 1:d){
    grid[,i] <- rep(Axis, each = grid.size^(i-1))
  }

  return(grid)
}
