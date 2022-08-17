
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mypac

<!-- badges: start -->
<!-- badges: end -->

This package implement the estimator for the incidence cumulative
function that introduced in the thesis “Nonparametric
cumulative-incidence estimation with delayed entry”. This work focuses
on focuses on the illness- death semi competing-risk model, a three
state stochastic model that is commonly used to describe the progress of
individuals from initial state to a terminal state, such as death. In
this work, a novel non-parametric estimator of the cumulative incidence
function under right-censoring, delayed entry, and semi-competing risks
is suggested. The new estimator uses the conditional survival function
of age at diagnosis, given age at death. A local linear estimation
approach has be used for estimating this conditional survival function.
It is assumed that given age at death, the recruitment age has no
additional predictive power for age at diagnosis. This assumption allows
the use of both prevalent and incident cases which is our main objective
and novel contribution, given existing works. An extensive simulation
study was performed for demonstrating the finite-sample properties of
the proposed estimation procedures.

## Installation

You can install the development version of mypac from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ShovalMarton/mypac")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(mypac)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
n = 5000
dat_A111 <- mypac::dat_A111 %>% subset(T2 > R)
dat_n <- dplyr::sample_n(tbl = dat_A111,size = n,replace = TRUE)

resul =  mypac::CIF(T1 = dat_n$T1,T2 = dat_n$T2,C = dat_n$C,R= dat_n$R,S_LT=0.97,Boot = 3,BW_method = "RT")
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
