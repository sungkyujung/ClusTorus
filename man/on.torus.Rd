% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/on.torus.R
\name{on.torus}
\alias{on.torus}
\title{Transform the angular values to be on principal interval}
\usage{
on.torus(x)
}
\arguments{
\item{x}{Matrix-formed numeric data which has 2 columns.}
}
\description{
\code{on.torus()} transforms data to be on \eqn{[0, 2\pi)^2}.
}
\examples{
\dontrun{
x <- matrix(c(10/3 * pi, 5/4 * pi), ncol = 2, byrow = TRUE)

on.torus(x)
}
}