#' Estimation of the Cumulative Incidence Function, CIF by two methods: KM and MG.
#' that represented on ???
#' Given T1, T2, C, R.
#' `creating_base_data()` returns database that contain that needed to the analysis. 
#' @param T2 The vector that represent the death age for each participant.
#'   The length of `T2` should be as the length of `T1`,`C`, `R`.
#'   Expected to be numeric vector.
#' @param T1 The vector that represent the disease age for each participant.
#'   The length of `T1` should be as the length of `T2`,`C`, `R`.
#'   Expected to be numeric vector.
#' @param C The vector that represent the censor age for each participant.
#'   The length of `C` should be as the length of `T2`,`T1`, `R`.
#'   Expected to be numeric vector.
#' @param R The vector that represent the recruited age for each participant.
#'   The length of `C` should be as the length of `T2`,`T1`, `C`.
#' @return data frame that contain:  
#'  T1,
#'  T2,
#'  C,
#'  V1 = pmin({{T1}},{{T2}},{{C}}),
#'  V2 = pmin(T2,C),
#'  delta1 = ifelse(V1==T1 ,1,0),
#'  delta2 = ifelse(V1==T2,1,0),
#'  delta3 = ifelse(T2<=C,1,0)*delta1,
#'  delta4 = ifelse(T2<=C,1,0)
#' @export
creating_data_base <- function(T1, T2, C, R,digits=0){

  # Check the input ----
  T1 <- as.numeric(T1)
  T2 <- as.numeric(T2)
  C <- as.numeric(C)
  R <- as.numeric(R)


  main_data <- data.frame(cbind(T1,T2,C,R))

  main_data <- dplyr::mutate(main_data,
                             V1 = pmin({{T1}},{{T2}},{{C}}),
                             V2 = pmin(T2,C),
                             delta1 = ifelse(V1==T1 ,1,0),
                             delta_1 = ifelse(T1 < T2,1,0),
                             delta2 = ifelse(V1==T2,1,0),
                             delta3 = ifelse(T2 < C,1,0)*delta1,
                             delta4 = ifelse(T2 < C,1,0),
                             delta_strike = as.numeric(delta1==1|delta4==1)
  )
  
  if(any(T2<R)){

    not_logical_val <- sum(T2<R)
    warning(paste('The following statment is satisfy: T2<R for',not_logical_val, 'rows'))
    #main_data <- subset(main_data,!(T2<R|C<R|T1>T2))
  }
  #main_data <- subset(main_data,T2>=R)
  main_data <- round(main_data,digits = digits)
  return(main_data)
}
