# ClusTorus

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ClusTorus)](https://cran.r-project.org/package=ClusTorus)
[![R-CMD-check](https://github.com/sungkyujung/ClusTorus/workflows/R-CMD-check/badge.svg)](https://github.com/sungkyujung/ClusTorus/actions)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![downloads](https://cranlogs.r-pkg.org/badges/ClusTorus)](https://cranlogs.r-pkg.org/badges/ClusTorus)
<!-- badges: end -->
Prediction and Clustering on The Torus by Conformal Prediction

This R package provides R functions to perform clustering multivariate angular data of any dimension on general tori. 

A part of the functions and data sets in this package was used to develop an unpublished paper "Clustering on the torus by conformal prediction" by Jung, Park and Kim. 

## Installation

```r
# install from CRAN
install.packages("ClusTorus")

# development version is available via this code:
# install.packages("devtools")
devtools::install_github("sungkyujung/ClusTorus")
```
If you also want to install the vignettes, then use the code:
```r
devtools::install_github("sungkyujung/ClusTorus", build_vignettes = TRUE)
```
