% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMsinvMmix.R
\name{EMsinvMmix}
\alias{EMsinvMmix}
\title{Fitting mixtures of bivariate von Mises distribution}
\usage{
EMsinvMmix(
  data,
  J = 4,
  parammat = EMsinvMmix.init(data, J),
  THRESHOLD = 1e-10,
  maxiter = 100,
  type = c("circular", "axis-aligned", "general"),
  kmax = 500,
  verbose = TRUE
)
}
\arguments{
\item{data}{n x 2 matrix of toroidal data on \eqn{[0, 2\pi)^2}}

\item{J}{number of components of mixture density}

\item{parammat}{6 x J parameter data with the following components:

  \code{parammat[1, ]} : the weights for each von Mises sine density

  \code{parammat[n + 1, ]} : \eqn{\kappa_n} for each von Mises
  sine density for n = 1, 2, 3

  \code{parammat[m + 4, ]} : \eqn{\mu_m} for each von Mises
   sine density for m = 1, 2}

\item{THRESHOLD}{number of threshold for difference between updating and
updated parameters.}

\item{maxiter}{the maximal number of iteration.}

\item{type}{a string one of "circular", "axis-aligned", "general",
and "Bayesian" which determines the fitting method.}

\item{kmax}{the maximal number of kappa. If estimated kappa is
larger than \code{kmax}, then put kappa as \code{kmax}.}

\item{verbose}{boolean index, which indicates whether display
additional details as to what the algorithm is doing or
how many loops are done.}
}
\value{
returns approximated parameters for bivariate normal
distribution with \code{list}:

  \code{list$Sigmainv[j]} : approximated covariance matrix for
  j-th bivariate normal distribution, approximation of the j-th von Mises.

  \code{list$c[j]} : approximated \eqn{|2\pi\Sigma|^{-1}} for
  j-th bivariate normal distribution, approximation of the j-th von Mises.
}
\description{
\code{EMsinvMmix} returns fitted parameters of J-mixture of
  bivariate sine von Mises distributions.
}
\details{
This algorithm is based on ECME algorithm. That is,
  constructed with E - step and M - step and M - step
  maximizes the parameters with given \code{type}.

  If \code{type == "circular"}, then the mixture density is
  just a product of two independent von Mises.

  If \code{type == "axis-aligned"}, then the mixture density is
  the special case of \code{type == "circular"}: only need to
  take care of the common concentration parameter.

  If\code{type == "general"}, then the fitting the mixture
  density is more complicated than before, check the detail of
  the reference article.
}
\examples{
\donttest{
data <- ILE[1:200, 1:2]

EMsinvMmix(data, J = 3,
           THRESHOLD = 1e-10, maxiter = 200,
           type = "general", kmax = 500, verbose = FALSE)
}
}
\references{
Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
}
