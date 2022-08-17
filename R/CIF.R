#' Estimation of the Cumulative Incidence Function. Gorfine's method,
#' as discribed by... and A 
#' Given T1, T2, C, R
#' `CIF()` return the estimated CIF from the imported data.
#'
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
#' @param kernel The kernel function that will be used in "PE" method.
#' Expected to be one of the follow: CosK(x) EpaK(x) Epa2K(x) gaussK(x), from "locpol" package see there for more details.
#' the default set to gaussK(x).
#' @param Boot The number of iterations in the bootstrap method. Default set to 100
#' @param S_LT The survival to live until the left truncation age, LT.
#'  The default calculated base on the data that the UK government published.
#'   @seealso https://en.wikipedia.org/wiki/Kernel_(statistics)#cite_note-4
#' @param digits specified number of decimal places (default 0) of `T1`,`C`, `R` and `T2`. 
#' @param LT the left truncation value, default set to 40
#' @param BW_method Bandwidth selection method, expected to be one of the follow: RT -  for rule of thumb,  plugin,
#'  for plugin method and  CV for cross validation. The recommended method is RT. These options are based on functions from "locpol" package.
#' @return list object that contains three objects: 'AJ_CIF' contains the CIF that estimated according to Aalen-Johansen method.
#'   'PE_CIF' contains the CIF that estimated according to MG method.
#'   'timestocheck' contains the responds ages. 
#' @export 
CIF = function(T1,T2,C,R,S_LT=NULL,BW_method="RT",Boot = 100,kernel = locpol::gaussK,digits = 0 ,LT=40){
  
  if(is.null(S_LT)){S_LT = (1-sum(death_group$p_t2[1:(LT+1)]))}
  
  datt = data.frame(cbind(T1,T2,C,R))
  
  datt <- mypac::creating_data_base(T1 = datt$T1,
                                    T2 = datt$T2,
                                    C = datt$C,
                                    R = datt$R,
                                    digits = 0)
  
  timestocheck <- 1:max(round(datt$T1[datt$delta1==1]))
  
  n = nrow(datt)
  
  Nelson.results = NULL
  My.results_all_External3 = NULL
  
  if(Boot>0){
    
    for(iter in 1:Boot){
      
      main_data_boot <- dplyr::sample_n(tbl = datt,size = n,replace = TRUE) 
      
      AJ_fun <- CIF_Aalen_Johansen(main_data = main_data_boot,digits = 0)
      
      
      new_fun <- CIF_new_kernel(main_data = main_data_boot,S_LT = S_LT,digits = 0,LT = LT,BW_method = BW_method)
      
      
      Nelson.results <- cbind(Nelson.results,AJ_fun(timestocheck))
      My.results_all_External3 <- cbind(My.results_all_External3,new_fun(timestocheck))
    }
    
    All.results <- list()
    
    All.results$AJ_CIF = Nelson.results
    All.results$PE_CIF = My.results_all_External3
    All.results$timestocheck = timestocheck
    
    return(All.results)
    
    
  }
  
  
  if(Boot == 0){
    
    AJ_fun <- CIF_Aalen_kernel(main_data = datt,digits = 0)
    new_fun <- CIF_new_kernel1(main_data = datt,S_LT = S_LT,digits = 0,LT = 40)

    
    Nelson.results <- cbind(Nelson.results,AJ_fun(timestocheck))
    My.results_all_External3 <- cbind(My.results_all_External3,new_fun(timestocheck))
    
    All.results <- list()
    
    All.results$AJ_CIF = Nelson.results
    All.results$PE_CIF = My.results_all_External3
    return(All.results)
  }
  
  
  
  
}

