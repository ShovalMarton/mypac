#' Simulated data set, for  T1=A, T2=1, C=1 and R=1.
#'
#' A simulated dataset containing  5*10^6 rows, that used in the thesis "Nonparametric cumulative-incidence
#' estimation with delayed entry"
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{T1}{Age at diagnosis, T1 ~ TWeibull(lambda=4,k=115,a=40)}
#'   \item{T2}{Age at death. T2 was generated from the empirical distribution of death age in the United Kingdom.
#'    For observations that were diagnosed with the disease, T1 < T2,
#'    the original $T_2$ value was discarded, and given T1, a new T2 value was sampled from a truncated Weibull distribution.
#'    For more details, see in the thesis "Nonparametric cumulative-incidence estimation with delayed entry"}
#'   \item{R}{Age at recruitment, R ~  U(40,69)}
#'   \item{C}{Age at censoring, C ~ U(0.1,35)}
#'  
#' }
"dat_A111"

#' Simulated data set, for T1=B, T2=1, C=1 and R=1.
#'
#' A simulated dataset containing  5*10^6 rows, that used in the thesis "Nonparametric cumulative-incidence
#' estimation with delayed entry"
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{T1}{Age at diagnosis, T1 ~ TWeibull(lambda=4,k=130,a=40)}
#'   \item{T2}{Age at death. T2 was generated from the empirical distribution of death age in the United Kingdom.
#'    For observations that were diagnosed with the disease, T1 < T2,
#'    the original $T_2$ value was discarded, and given T1, a new T2 value was sampled from a truncated Weibull distribution.
#'    For more details, see in the thesis "Nonparametric cumulative-incidence estimation with delayed entry"}
#'   \item{R}{Age at recruitment, R ~  U(40,69)}
#'   \item{C}{Age at censoring, C ~ U(0.1,35)}
#'  
#' }
"dat_B111"



#' Simulated data set, for T1=C, T2=1 C=1 and R=1.
#'
#' A simulated dataset containing  5*10^6 rows, that used in the thesis "Nonparametric cumulative-incidence
#' estimation with delayed entry"
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{T1}{Age at diagnosis, T1 ~ Weibull(lambda=3.5,k=20)}
#'   \item{T2}{Age at death. T2 was generated from the empirical distribution of death age in the United Kingdom.
#'    For observations that were diagnosed with the disease, T1 < T2,
#'    the original $T_2$ value was discarded, and given T1, a new T2 value was sampled from a truncated Weibull distribution.
#'    For more details, see in the thesis "Nonparametric cumulative-incidence estimation with delayed entry"}
#'   \item{R}{Age at recruitment, R ~  U(40,69)}
#'   \item{C}{Age at censoring, C ~ U(0.1,35)}
#'  
#' }
"dat_C111"


###########
#' Simulated data set, for T1=A, T2=1 C=2 and R=1.
#'
#' A simulated dataset containing  5*10^6 rows, that used in the thesis "Nonparametric cumulative-incidence
#' estimation with delayed entry"
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{T1}{Age at diagnosis, T1 ~ TWeibull(lambda=4,k=115,a=40)}
#'   \item{T2}{Age at death. T2 was generated from the empirical distribution of death age in the United Kingdom.
#'    For observations that were diagnosed with the disease, T1 < T2,
#'    the original $T_2$ value was discarded, and given T1, a new T2 value was sampled from a truncated Weibull distribution.
#'    For more details, see in the thesis "Nonparametric cumulative-incidence estimation with delayed entry"}
#'   \item{R}{Age at recruitment, R generated from the empirical distribution of the UKB recruitment age.}
#'   \item{C}{Age at censoring, C ~ U(11,15)}
#'  
#' }
"dat_A122"

#' Simulated data set, for T1=B, T2=1, C=2 and R=2.
#'
#' A simulated dataset containing  5*10^6 rows, that used in the thesis "Nonparametric cumulative-incidence
#' estimation with delayed entry"
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{T1}{Age at diagnosis, T1 ~ TWeibull(lambda=4,k=130,a=40)}
#'   \item{T2}{Age at death. T2 was generated from the empirical distribution of death age in the United Kingdom.
#'    For observations that were diagnosed with the disease, T1 < T2,
#'    the original $T_2$ value was discarded, and given T1, a new T2 value was sampled from a truncated Weibull distribution.
#'    For more details, see in the thesis "Nonparametric cumulative-incidence estimation with delayed entry"}
#'   \item{R}{Age at recruitment, R generated from the empirical distribution of the UKB recruitment age.}
#'   \item{C}{Age at censoring, C ~ U(11,15)}
#'  
#' }
"dat_B122"



#' Simulated data set, for T1=C, T2=1, C=2 and R=2, .
#'
#' A simulated dataset containing  5*10^6 rows, that used in the thesis "Nonparametric cumulative-incidence
#' estimation with delayed entry"
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{T1}{Age at diagnosis, T1 ~ Weibull(lambda=3.5,k=20)}
#'   \item{T2}{Age at death. T2 was generated from the empirical distribution of death age in the United Kingdom.
#'    For observations that were diagnosed with the disease, T1 < T2,
#'    the original $T_2$ value was discarded, and given T1, a new T2 value was sampled from a truncated Weibull distribution.
#'    For more details, see in the thesis "Nonparametric cumulative-incidence estimation with delayed entry"}
#'   \item{R}{Age at recruitment, R generated from the empirical distribution of the UKB recruitment age.}
#'   \item{C}{Age at censoring, C ~ U(11,15)}
#'  
#' }
"dat_C122"
