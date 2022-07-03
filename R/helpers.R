
#' @export
CIF_new_kernel <- function(main_data,kernel = locpol::gaussK,digits = 0, 
                           S_LT, LT=LT,BW_method="RT"){
  
  # Check the input ----
  T2_observed <- sort(unique(round(main_data$V2[main_data$delta4==1])))
  T1_observed <- sort(unique(round(main_data$T1[main_data$delta1==1])))
  
  main_data <- round(main_data)
  
  # Orginaized data ----
  main_data <- dplyr::mutate(main_data, death_time = V1*delta2+V2*delta3)#,death_indi=delta2+delta3)
  
  main_data$death_time <- apply(main_data,1,function(x){ifelse(x["delta2"]==1 & x["delta3"]==1,
                                                               x["V1"]*x["delta2"],x["death_time"])})
  
  T2_observed_n <-  length(T2_observed)
  T1_observed_n <- length(T1_observed)
  
  main_data_b_sub <- subset(main_data, death_time > 0)
  
  # 1: creating the hazard function  ----
  if(BW_method == "RT"){
    est.int1 <-  H_fun_BW_RT(main_data_b_sub = main_data_b_sub, T1_observed = T1_observed,T2_observed = T2_observed,kernel_fun = kernel)
  }
  if(BW_method == "RT_1000"){
    est.int1 <-  H_fun_BW_RT_1000(main_data_b_sub = main_data_b_sub, T1_observed = T1_observed,T2_observed = T2_observed,kernel_fun = kernel)
  }
  if(BW_method == "RT_01"){
    est.int1 <-  H_fun_BW_RT_01(main_data_b_sub = main_data_b_sub, T1_observed = T1_observed,T2_observed = T2_observed,kernel_fun = kernel)
  }
  if(BW_method == "plugin"){
    est.int1 <-  H_fun_BW_plugin(main_data_b_sub = main_data_b_sub, T1_observed = T1_observed,T2_observed = T2_observed,kernel_fun = kernel)
  } 
 
  if(BW_method == "CV"){
    est.int1 <-  H_fun_BW_CV(main_data_b_sub = main_data_b_sub, T1_observed = T1_observed,T2_observed = T2_observed,kernel_fun = kernel)
  } 
  
  # 2: creating the given CH  ----
  
  Lambda.given.s <- matrix(NA,T2_observed_n,T1_observed_n)
  
  for (ii in 1:T2_observed_n) { #s
    for (jj in 1:T1_observed_n) { #t1, t1<=s
      
      Lambda.given.s[ii,jj] <- sum(est.int1[ii,1:jj],na.rm = TRUE)
    }
  }
  Surv.given.s <- exp(-(Lambda.given.s))
  # 3: Creating the survival function ----
  
  # Estimation of T2 using KM approach, the estimation is entirely from the observe data.
  main_data[which(main_data$V2==main_data$R),"V2"]=main_data[which(main_data$V2==main_data$R),"V2"]+0.5
  
  fit.3a <- summary(survival::survfit(survival::Surv(time  = main_data$R,time2 =  main_data$V2,
                                                     event = main_data$delta4,type = "counting")~1))
 
  time.3a <- fit.3a$time; h.3a <- fit.3a$n.event/fit.3a$n.risk#; ch.3a <- cumsum(h.3a);
  S.3a <-  cumprod(c(1-h.3a))
  F.3a <- stepfun(time.3a,c(0,1-S.3a))

  val_all <- F.3a(c(min(T2_observed)-1,T2_observed))
  dval_all<-c(diff(val_all,lag=1))
  
  # 4: Estimation of T2 using external data.
  
  dval_half_Malka <- dval_all*S_LT
  
  # 5: Final estimator.
  Final.est <- c()
  
  for (jj in 1:T1_observed_n) {
    Final.est[jj] <- sum((1-Surv.given.s[,jj])*dval_half_Malka,na.rm = TRUE)
  }
  
  Final.est  = stepfun(T1_observed,c(0,Final.est))
  
  return(Final.est)
  
}

#' @export
CIF_new_kernel_Cpp <- function(main_data,kernel = 1, S_LT = NULL){#},LT=LT){

  # Check the input ----
  main_data_b_sub <- subset(main_data, death_time > 0 & death_time != Inf & R <= V2) %>% as.matrix()
  
  T2_observed <- sort(unique(round(main_data$V2[main_data$delta4==1])))
  T1_observed <- sort(unique(round(main_data$T1[main_data$delta1==1])))
  #t1 = s
  
  est.int2 =  mypac:::est_h(data  = main_data_b_sub,krnl = 1,S_s = T2_observed ,T1_s = T1_observed)


  # 2: creating the given survival  ----

  Rcpp_cum_haz_given =  mypac:::Lambda_given_s(est_int1 = est.int2,S_s = T2_observed ,T1_s = T1_observed)
  Rcpp_surv_given = exp(-Rcpp_cum_haz_given)
  
  # 3: Creating the survival function ----
  # Estimation of T2 using KM approach, the estimation is entirely from the observe data.
  
  
  main_data[which(main_data$V2==main_data$R),"V2"] = main_data[which(main_data$V2==main_data$R),"V2"]+0.5
  
  fit.3a <- summary(survival::survfit(survival::Surv(time  = main_data$R,time2 =  main_data$V2, event = main_data$delta4,type = "counting")~1))
  time.3a <- fit.3a$time; h.3a <- fit.3a$n.event/fit.3a$n.risk;
  
  S.3a <-  cumprod(c(1-h.3a))
  F.3a <- stepfun(time.3a,c(0,1-S.3a))
 
  val_all <- F.3a(c(min(T2_observed)-1,T2_observed))
  dval_all<-c(diff(val_all,lag=1))
  

  # 4: Estimation of T2 using external data. ----

  dval_half_Malka <- dval_all*S_LT

  # 5: Final estimator.-------

  Final.est =  mypac:::final_CIF(Surv_given_s = Rcpp_surv_given ,T1_s = T1_observed,dval = dval_half_Malka)
  estG_tilde <- stepfun(T1_observed,c(0,Final.est))
  return(estG_tilde)#Final.est)

}

#' @export
CIF_Aalen_Johansen = function(main_data,digits = 0){
  
  # Orgnazie the data----
  
  NAE_data =  round(subset(main_data,V1>R))
  NAE_data[which(NAE_data$V1==NAE_data$R),"V1"] = NAE_data$V1[which(NAE_data$V1==NAE_data$R)]+0.5
  
  
  fit.12 <- summary(survival::survfit(survival::Surv(time = NAE_data$R, time2 = NAE_data$V1,
                                                     event = NAE_data$delta1,type = "counting")~1))
  time.12 <- fit.12$time; h.12 <- fit.12$n.event/fit.12$n.risk; ch.12 <- cumsum(h.12);
  chstepfun.12 <- stepfun(time.12,c(0,ch.12))
  
  
  #2---------
  fit.strike <- summary(survival::survfit(survival::Surv(time = NAE_data$R, time2 = NAE_data$V1,
                                                         event = NAE_data$delta_strike,type = "counting")~1))
  time.strike<- fit.strike$time
  
  S_stepfun.strike <- stepfun(time.strike,c(1,fit.strike$surv))
  S.time.strikeat12 <- S_stepfun.strike(time.12)
  n_12 <- length(time.12)
  
  estG <- cumsum(h.12*c(1,S.time.strikeat12[-n_12]))
  estG <- stepfun(time.12,c(0,estG))
  
  return(estG)
  
}

#' @export
groups_division <- function(x,divide,groups){
  n_groups = 70/divide - 40/divide
  groups =  40+divide*c(1:n_groups)
  return(groups[which(x["R"]<=groups)[1]])
}

#' @export
H_fun_BW_RT <- function(main_data_b_sub, T1_observed,T2_observed,kernel_fun){
 
  
  T2_observed_n  <-  length(T2_observed)
  T1_observed_n <- length(T1_observed)
  
  est.int1 <- matrix(NA,T2_observed_n,T1_observed_n)
  #mat_bw <- matrix(NA,T2_observed_n,T1_observed_n)
  bw1=NA
  
  for (ii in 1:T2_observed_n) { #s
    for (jj in 1:T1_observed_n) { #t1, t1<=s
      
      if(T2_observed[ii] < T1_observed[jj]){next()}
      
      Yvec.temp <- ifelse(main_data_b_sub[,"V1"]>=T1_observed[jj],1,0)#
      Nvec <- ifelse((main_data_b_sub[,"V1"]==T1_observed[jj])&(main_data_b_sub[,"delta1"]==1),1,0)
      Yvec <- subset(Yvec.temp,Yvec.temp>0)
      Nvec <- subset(Nvec,Yvec.temp>0)
      y.var <- Nvec/Yvec
      x.var <- subset(main_data_b_sub[,"death_time"],Yvec.temp>0) - T2_observed[ii]
      
      #RT bandwidth locpol
      #TODO write a a function by myself, it will be easier with  Rcpp.
      try(bw1 <- locpol::thumbBw(x = x.var, y = y.var, deg = 1, kernel = kernel_fun,
                                 weig = rep(1,length(y.var))))
      
      
      if(any(is.nan(bw1)|is.na(bw1)|length(bw1)==0)){next() }
      
      weight1 <- kernel_fun(x = x.var/bw1)
      #mat_bw[ii,jj] <- bw1
      try(fit.lm1 <- lm(y.var ~ x.var, weights = weight1))
      est.int1[ii,jj] <- fit.lm1$coefficients[1]
      
      bw1=NA
    }
  }
  return(est.int1)
  
}

#' @export
H_fun_BW_CV <- function(main_data_b_sub, T1_observed,T2_observed,kernel_fun){
  
  
  T2_observed_n <-  length(T2_observed)
  T1_observed_n <- length(T1_observed)
  
  est.int1 <- matrix(NA,T2_observed_n,T1_observed_n)
 
  bw1=NA
  
  tictoc::tic()
  for (ii in 1:T2_observed_n) { #s
    for (jj in 1:T1_observed_n) { #t1, t1<=s
      
      if(T2_observed[ii]<T1_observed[jj]){next()}
      
      Yvec.temp <- ifelse(main_data_b_sub[,"V1"]>=T1_observed[jj],1,0)#
      Nvec <- ifelse((main_data_b_sub[,"V1"]==T1_observed[jj])&(main_data_b_sub[,"delta1"]==1),1,0)
      Yvec <- subset(Yvec.temp,Yvec.temp>0)
      Nvec <- subset(Nvec,Yvec.temp>0)
      y.var <- Nvec/Yvec
      x.var <- subset(main_data_b_sub[,"death_time"],Yvec.temp>0) - T2_observed[ii]
      
      #RT bandwidth locpol
      #TODO write a a function by myselef, it will be easyer with  Rcpp.
      try(bw1 <- locpol::regCVBwSelC(x = x.var, y = y.var, deg = 1, kernel = kernel_fun,
                                 weig = rep(1,length(y.var))))
      
      
      if(any(is.nan(bw1)|is.na(bw1)|length(bw1)==0)){next() }
    
      weight1 <- kernel_fun(x.var/bw1)
      try(fit.lm1 <- lm(y.var ~ x.var, weights = weight1))
      est.int1[ii,jj] <- fit.lm1$coefficients[1]
      
      bw1=NA
    }
  }
  tictoc::toc()
  
  return(est.int1)
  
}

#' @export
H_fun_BW_plugin <- function(main_data_b_sub, T1_observed,T2_observed,kernel_fun){
  
  
  T2_observed_n <-  length(T2_observed)
  T1_observed_n <- length(T1_observed)
  
  est.int1 <- matrix(NA,T2_observed_n,T1_observed_n)
  bw1=NA
  
  for (ii in 1:T2_observed_n) { #s
    for (jj in 1:T1_observed_n) { #t1, t1<=s
      
      if(T2_observed[ii]<T1_observed[jj]){next()}
      
      Yvec.temp <- ifelse(main_data_b_sub[,"V1"]>=T1_observed[jj],1,0)#
      Nvec <- ifelse((main_data_b_sub[,"V1"]==T1_observed[jj])&(main_data_b_sub[,"delta1"]==1),1,0)
      Yvec <- subset(Yvec.temp,Yvec.temp>0)
      Nvec <- subset(Nvec,Yvec.temp>0)
      y.var <- Nvec/Yvec
      x.var <- subset(main_data_b_sub[,"death_time"],Yvec.temp>0) - T2_observed[ii]
      
      #RT bandwidth locpol
      #TODO write a a function by myselef, it will be easyer with  Rcpp.
      try(bw1 <- locpol::pluginBw(x = x.var, y = y.var, deg = 1, kernel = kernel_fun,
                                     weig = rep(1,length(y.var))))
      
      
      if(any(is.nan(bw1)|is.na(bw1)|length(bw1)==0)){next() }
      
      weight1 <- kernel_fun(x.var/bw1)
      try(fit.lm1 <- lm(y.var ~ x.var, weights = weight1))
      est.int1[ii,jj] <- fit.lm1$coefficients[1]
      
      bw1=NA
    }
  }
  return(est.int1)
  
}


#' `CIF_new_assumption_checking()` returns assumption test graph, as it describe in ....
#' @param T2 The vector that represents the death age for each participant.
#'   The length of `T2` should be as the length of `T1`,`C`, `R`.
#'   Expected to be numeric vector.
#'
#' @param T1 The vector that represents the disease age for each participant.
#'   The length of `T1` should be as the length of `T2`,`C`, `R`.
#'   Expected to be numeric vector.
#' @param C The vector that represents the censor age for each participant.
#'   The length of `C` should be as the length of `T2`,`T1`, `R`.
#'   Expected to be numeric vector.
#' @param R The vector that represents the recruited age for each participant.
#'   The length of `R` should be as the length of `T2`,`T1`, `C`.
#'   Expected to be numeric vector.
#' @param kernel The kernel function that will be used in "MG" method.
#' Expected to be one of the follow: 1 - Gaussian, 2 - Epanechnikov,3 - box, 4 - Logistic, 5 - Sigmoid function.
#' the default set to 1.@seealso https://en.wikipedia.org/wiki/Kernel_(statistics)#cite_note-4
#' @param S_LT The survival to live until the left truncation age, LT.
#'  The default calculated base on the data that the UK government published.
#' @param LT The left truncation data of the imported data. If NULL the Left truncation will be the minimum recruited
#'  age in the risk set.
#' @param division_val The gap of years in it the estimator will be examine. Foe example, for division_val=5,
#'   when the range of recruitment is 40-70, every 5 years the estimator performance will be shown.
#' @return a list that contains the objects: ggplot object with graph, that contains multiple curves, one for each subgroups;
#'  Martix that contains the AIC for each object and for whole data.
#' @export
CIF_new_assumption_checking <- function(T1,T2,C,R,S_LT, division_val = NULL, kernel = locpol::gaussK, LT=40){
  main_data <- creating_data_base(T1,T2,C,R)
  
  if(is.null(division_val)){division_val = 5}
  n_groups = 70/division_val - LT/division_val
  main_data$group <- apply(main_data,1,mypac:::groups_division,divide = division_val)
  main_data_List <- split(main_data,as.factor(main_data$group))
  
  final_data_est <- NULL
  final_label <- c()
  
  for(i in 1:length(main_data_List)){
    
    main_data_i <- main_data_List[[i]]
    timestocheck <- unique(round(main_data_i$T1[main_data_i$delta1==1]))
    
    results = CIF_new_kernel(main_data = main_data_i,S_LT = S_LT,LT = LT)
    Final.est_i = results(timestocheck)
    
    data_final_i <- data.frame(cbind(Age_T1 = as.numeric(timestocheck),
                                     Final.est_i = as.numeric(round(Final.est_i,10)),
                                     label = max(main_data_i$R)))
    
    final_label[i] <- paste("R = " ,min(main_data_i$R),"-",max(main_data_i$R),sep = "")
    final_data_est <-  rbind(final_data_est,data_final_i)
    
  }
  
  
  xlab <- expression(paste(Age))
  
  final_data_est <- data.frame(final_data_est)
  final_data_est <- final_data_est %>% arrange(label,Age_T1)
  final_data_est$label <- factor(final_data_est$label,levels = unique(final_data_est$label),labels = final_label)
 
  q <- ggplot2::ggplot(final_data_est, aes(x = Age_T1, y = Final.est_i)) +
    geom_step(size=1.5) +
    facet_grid(~ label)+
    labs(y= "Final estimaor", x = xlab, sep="") +
    #ggtitle(label = "Checking S(T1|T2) = S(T1|T2,R)",) +
    theme(axis.title = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 20, face = "bold"), 
          strip.text.x = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12))+
    scale_x_continuous(sec.axis = sec_axis(~ .*10 , breaks = NULL, labels = NULL))#, name = "Recruitment Age"
  
  list_diagnosis <-  list()
  list_diagnosis$Plot_Diagnosis <- q
  
  return(list_diagnosis)
}