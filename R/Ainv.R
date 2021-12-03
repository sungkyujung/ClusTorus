# Inverse of \eqn{A}
#
# \code{Ainv} returns the inverse of A, where A is the ratio of the first
#   and zeroth order Bessel functions of the first kind with finite series
#   approximation. If the inverse value \eqn{k} of \eqn{A(x)} is larger than
#   \code{kmax}, then it returns \code{kmax}.
#
# @param x numeric value in \eqn{[0,1]}.
# @param kmax positive value. Default value is 2000.
# @seealso See the similar function in package \code{circular} :
#   \code{\link[circular]{A1inv}}
# @references  Mardia, K. V. and Jupp, P. E. (2000), \emph{Directional Statistics},
#   New York: Wiley.
# @examples
# \dontrun{
# Ainv(0.7, kmax = 3000)
# }
Ainv <- function (x, kmax = 2000){
  if (!is.numeric(x) || length(x) > 1) {stop("Invalid value to x : x must be in the unit interval.")}
  if (!is.numeric(kmax) || length(kmax) > 1) {stop("Invalid value to kmax : kmax must be a numeric scalar.")}
  k <- ifelse(0 <= x & x < 0.53, 2 * x + x^3 + (5 * x^5)/6,
              ifelse(x < 0.85, -0.4 + 1.39 * x + 0.43/(1 - x),
                     ifelse(x < 1,1/(x^3 - 4 * x^2 + 3 * x), kmax)))
  ifelse(k > kmax, kmax, k)
}
