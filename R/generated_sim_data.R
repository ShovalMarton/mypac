#' Generate simulated data for the death-illness model
#' `sim_data_fun()` returns data frame that contains T1,T2,C and R.
#' To let the user explore the estimators' behavior in
#'  different types of diseases, recruitment processes, and censoring.
#'
#' @param n number of observations. The length(n) expected to be 1.
#'
#' @param R_method The type of distribution that the vector R will generate from.
#'  expected to be one of the numbers 0-2.
#'
#'  Type 0: set all the observations to be R=LT.
#'  Type 1: generate the observations to from discrete uniform distribution, between LT to U_R.
#'  Type 2: generate the observations from the UKB's recruitment distribution.
#'
#' @param death_method The death distribution is based on the UK population.
#'  death_method parameter refers to the disease fatality. Expected to values between 0-3.
#'  Type 0 - The disease does not affect life expectancy.
#'  Types 1-3 represent different disease fatalities, with dependence on the disease age.
#'  All these types are based on the truncated Weibull distribution (Tweibull).
#'  1. Type 1: Tweibull(a = T1,b = 105,shape = 5.5 ,scale=55)
#'  2. Type 2: Tweibull(a = T1,b = 105,shape = 5.5 ,scale=68)
#'  3. Type 3: Tweibull(a = T1,b = 105,shape = 5.5 ,scale= 75)
#'  
#' @param T1_method The disease distribution. Expected to values A-C.
#'  1. Type A: Tweibull(a = 40, shape = 4 ,scale=115)
#'  2. Type B: Tweibull(a = 40, shape = 4 ,scale=130)
#'  3. Type C: Weibull(a = T1, shape=3.5, scale=200)
#'
#' @param LT The left truncation age of the simulated study,
#'  default set to be 40 as the UKB data.
#'
#' @param U_R The upper edge that observation can get into the simulated study.
#'  default set to be 69 as the UKB data.
#' @return The function returns a data frame with simulated T1, T2, C, R.
#' @export
sim_data_fun <- function(n,R_method,death_method,T1_method,C_method,LT=40, U_R=69){

  if(length(n) != 1){stop('length(n) must be equal to 1')}
  if(length(R_method) != 1){stop('length(R) must be equal to 1')}
  if(!c(T1_method %in% c("A","B","C"))){stop('T1_method expected to be one of the following arguments: "A","B","C"')}
  if(!c(death_method %in% c(0:3))){stop('death_method expected to be one of the following arguments: of 0,1,2,3')}
  if(!c(R_method %in% c(0:2))){stop('R_method expected to one of the following arguments: 0,1,2')}

  # Generate T1 ----
  if(T1_method=="A"){T1<- truncdist::rtrunc(n, spec="weibull",shape = 4 ,scale=115,a = 40)} #5
  if(T1_method=="B"){T1<- truncdist::rtrunc(n, spec="weibull",shape = 4 ,scale=130,a = 40)} #6
  #if(T1_method==16){T1 <- truncdist::rtrunc(n, spec="weibull",shape = 0.9 ,scale=420,a = 40)} 
  if(T1_method=="C"){T1<- rweibull(n, shape=3.5, scale=200)} #18
  
  # Generate T2 ----
  T2 <- abs(sample(c(0:105),size = n,prob = death_group$p_t2,replace = TRUE) + runif(n,-0.5,0.5))

  data <- data.frame(cbind(T1,T2))
  data <- data %>%  dplyr::mutate(delta1 = T1 <= T2)

  data$T2 <- apply(X = data, MARGIN =1, FUN = mypac::update_T2, type=death_method)
 

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
  
  data <- data %>% dplyr::select(-delta1)

  data
}

#' @export
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





