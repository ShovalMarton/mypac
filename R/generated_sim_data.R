#' Generate simulated data for the death-illness model
#' `sim_data_fun()` returns data frame that contains T1,T2,C and R.
#'  In order to let the user explore the estimators' behavior in
#'  difference type of diseases and recruitment processes there are multiple options.
#'  The different type of data allow to explore how the fatality,
#'  ages' distribution of the eruption, death age and different dependency
#'  between the disease and the recruitment process affects on the estimators.
#'
#' @param n number of observations. The length(n) expected to be 1.
#'
#' @param R_method The type of distribution that the vector R will generate from.
#'  expected to be one of the numbers 0-7.
#'
#'  Type 0: set all the observations to be R=LT.
#'  Type 1: generate the observations to from discrete uniform distribution, between LT to U_R.
#'  Type 2: generate the observations from the UKB's recruitment distribution.
#'  Type 3-6: generate the observations from normal truncation distribution, between LT to U_R to be
#'  with dependency upon delta 1:
#'  1. Type 3
#'  if the observation i is sick (delta1==1), R_i will arrive from the distribution:
#'
#'  TN(a=LT, b=U_R, mean = 50, sd = 5),
#'
#'  otherwise from the distribution:
#'  TN(a=LT, b=U_R, mean = 55, sd = 10).
#'  2. Type 4
#'  if the observation i is sick, R_i will arrive from the distribution:
#'
#'  TN(a=LT, b=U_R, mean = 55, sd = 10),
#'
#'  otherwise from the distribution:
#'  TN(a=LT, b=U_R, mean = 50, sd = 5).
#'
#'  3. Type 5
#'  if the observation i is sick, R_i will arrive from the distribuation:
#'
#'  TN(a=LT, b=U_R, mean = 50, sd = 5),
#'
#'  otherwise from the distribution:
#'  TN(a=LT, b=U_R, mean = 60, sd = 7).
#'
#'
#'  4. Type 6
#'  if the observation i is sick, R_i will arrive from the distribuation:
#'
#'  TN(a=LT, b=U_R, mean = 60, sd = 7),
#'
#'  otherwise from the distribution:
#'  TN(a=LT, b=U_R, mean = 50, sd = 5).
#'
#'
#' @param death_method The death distrbution is base on the UK population.
#'  death_method parmater refers to the disease fatalty. Expect to recive values between 0-4.
#'  Type 0 - The disease not effect on the life expectancy.
#'  Types 1-4 are distrbutions  that represents diffrents disease fatality, with dependenct on the disease age.
#'  All this types base on the truncated normal distbution, when a set to be T1.
#'  1. Type 1: NT(a = T1,b = 105,sd=6,mean=70 )
#'  2. Type 2: NT(a = T1,b = 105,sd=7,mean=50 )
#'  3. Type 3: NT(a = T1,b = 105,sd=7,mean=20 )
#'  4. Type 4: NT(a = T1,b = 105,sd=6,mean=40 )
#'
#' @param T1_method Vector of the estimator's output, numeric vector.
#'  The length of CIF_new_esti expected to be equal to CIF_new_sd.
#'
#' @param LT The left truncation age of the simulate study,
#'  defult set to be 40 as the UKB data.
#'
#' @param U_R The upper adge that observation can get into the simulate study.
#'  defult set to be 69 as the UKB data.
#' @return The function return a data frame with simulate T1,T2,C,R.
#' @export
sim_data_fun <- function(n,R_method,death_method,T1_method,C_method,LT=40, U_R=69){

  if(length(n) != 1){stop('length(n) must be equal to 1')}
  #if(!c(T1_method %in% c(1:5))){stop('T1_method expected to be one of the follow arguments: 1,2,3,4,5')}
  #if(!c(death_method %in% c(0:4))){stop('death_method expected to be one of the follow arguments: of 0,1,2,3,4')}
  #if(!c(R_method %in% c(0:6))){warning('R_method expected to one of the follow arguments: 0,1,2,3,4,5,6 - R vector will not generte')}

  # Generate T1 ----
  if(T1_method==5){T1<- truncdist::rtrunc(n, spec="weibull",shape = 4 ,scale=115,a = 40)}
  if(T1_method==6){T1<- truncdist::rtrunc(n, spec="weibull",shape = 4 ,scale=130,a = 40)}
  if(T1_method==16){T1 <- truncdist::rtrunc(n, spec="weibull",shape = 0.9 ,scale=420,a = 40)}
  if(T1_method==18){T1<- rweibull(n, shape=3.5, scale=200)}
  
  # Generate T2 ----
  T2 <- abs(sample(c(0:105),size = n,prob = death_group$p_t2,replace = TRUE) + runif(n,-0.5,0.5))

  data <- data.frame(cbind(T1,T2))
  data <- data %>%  dplyr::mutate(delta1 = T1 <= T2)

  data$T2 <- apply(X = data, MARGIN =1, FUN = mypac::update_T2, type=death_method)
  #mean(data$T2[data$delta1==1] > data$T2_new1[data$delta1==1])
  
  #data$T2 <- ifelse(test = data$T2 > data$C,yes = Inf,no = data$T2)

  # Generate R ----
  if(R_method==0){ data$R = LT}
  if(R_method==1){ data$R = runif(n,LT,U_R)} #not depends on disease status
  if(R_method==2){ #not depends on disease status  - ukbb
    data$R =  sample(R_group$R,size = n,prob = R_group$p,replace = TRUE) + runif(n,-0.5,0.5)
  }

  # Generate C ----
  if(C_method == 1){data$C <- data$R + runif(n,0.1,35)}
  if(C_method == 2){data$C <- data$R + runif(n,11,15)}
  if(C_method == 3){data$C <- data$R + runif(n,0.1,15)}
  if(C_method == 4){data$C <-truncdist::rtrunc(n, spec="weibull",shape=9, scale=92,a = LT)}
  data <- data %>% dplyr::select(-delta1)

  data
}

update_T2 =  function(x,type=type){

  # In the older ages T1 occurs to increase in the expected age death.
  
  if(type==0){return( as.numeric(x["T2"]))}
  
  if(type==1){
    if(x["delta1"]==0|as.numeric(x["T2"])>=105){return( as.numeric(x["T2"]))}
    if(x["delta1"]==1){
      y = as.numeric(x["T1"])
      
      #return(min(truncdist::rtrunc(n = 1, spec="weibull",shape = 5.5 ,scale=55,a = y,b = 105),as.numeric(x["T2"])))
      return(min(truncdist::rtrunc(n = 1, spec="weibull",shape = 5.5 ,scale=55,a = y,b = 105)))
      
    }
    
  }
  
  
  if(type==2){
    if(x["delta1"]==0|as.numeric(x["T2"])>=105){return( as.numeric(x["T2"]))}
    if(x["delta1"]==1){
      y = as.numeric(x["T1"])

      return(min(truncdist::rtrunc(n = 1, spec="weibull",shape = 5.5 ,scale=68,a = y,b = 105)))#,as.numeric(x["T2"])

    }
  }


  if(type==3){
    if(x["delta1"]==0|as.numeric(x["T2"])>=105){return( as.numeric(x["T2"]))}
    if(x["delta1"]==1){
      y = as.numeric(x["T1"])

      return(min(truncdist::rtrunc(n = 1, spec="weibull",shape = 5.5 ,scale= 75,a = y,b = 105)))#,as.numeric(x["T2"])
      

    }
  }

}





