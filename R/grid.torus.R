#' Grid on torus
#'
#' \code{grid.torus()} returns an equally-spaced grid on torus.
#'
#' @param grid.size number of grid for each axis. Default value is 100.
#' @return \code{grid.torus} returns (grid.size) x (grid.size) numeric matrix
#'   which indicates the grid points on torus
#' @export
#' @examples
#' \dontrun{
#' grid.torus()
#' }
grid.torus <- function(grid.size = 100){
  # returns grid points on torus of size (grid.size) x (grid.size)
  Axis <- seq(-pi, pi, length = grid.size)
  return( cbind(rep(Axis, grid.size), rep(Axis, each=grid.size)))
}
