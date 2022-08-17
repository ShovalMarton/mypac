#' Export Graphs of CIF
#' `Graphs()` function returns plots of the CIF estimators.
#'  The output is a list og ggplot2 object, hence can be modified by adding ggplot objects using +
#'
#' @param years Vector of the corspoding years to the estimator,
#'   The length of `T2` should be as the length of `T1`,`C`, `R`.
#'   Expected to be numeric vector. If empty, the years will set to be
#'   years = c(1:length(CIF_new_sd)
#'
#' @param CIF_new_sd Vector of the estimator's standard deviation, numeric vector.
#'  The length of CIF_new_sd expected to be equal to CIF_new_esti
#'
#' @param CIF_new_esti Vector of the estimator's output, numeric vector.
#'  The length of CIF_new_esti expected to be equal to CIF_new_sd.
#' @export
Graphs =  function(years=NULL,CIF_new_sd,CIF_new_esti,CIF_nelson_sd=NULL,CIF_nelson_esti=NULL,alpha = 0.05){
if(is.null(years)){
  years = c(1:length(CIF_new_sd))
}
  plot_list <- list()
  #alpha <- 0.05

  if(!is.null(CIF_nelson_sd)){
    se.Nelson10<- CIF_nelson_sd
    mean.Nelson10 <- CIF_nelson_esti

    CI_NA_U <- mean.Nelson10+(1.96*(se.Nelson10))
    CI_NA_L <- mean.Nelson10-(1.96*(se.Nelson10))
  }

  mean.My10_all <- CIF_new_esti
  se.My10_all <- CIF_new_sd

  CI_NE_U_all <- mean.My10_all+(1.96*(se.My10_all))
  CI_NE_L_all <- mean.My10_all-(1.96*(se.My10_all))

  if(!is.null(CIF_nelson)){
    dat_graphs <- data.frame(years = CIF_new[,1],
                             mean.Nelson10,mean.My10_all,
                             se.Nelson10,se.My10_all,
                             CI_NA_U,CI_NA_L,
                             CI_NE_U_all,CI_NE_L_all)

    #--Define title and axis labels:
    main <- paste("The CIF,","n=",sample_size, sep = "")
    xlab <- expression(paste(Age))
    ylab <- expression(paste(G(Age)))

    #--Line styles and colours:
    col <- c("#1B9E77", "#D95F02")
    lty <- 1:2	# Line type
    lwd <- 3.5	# Line width


    plot_list$Naive_Proposed <- print(ggplot(dat_graphs,aes(x = years)) +
                                        geom_ribbon(aes(ymin=CI_NE_L_all,ymax=CI_NE_U_all),fill=col[1],alpha=0.3)+
                                        geom_ribbon(aes(ymin=CI_NA_L,ymax=CI_NA_U),fill=col[2],alpha=0.3)+

                                        geom_line(aes( y= mean.Nelson10,linetype ="Naive", col="Naive"), size=0.75)+
                                        geom_line(aes( y= mean.My10_all,linetype = "The Proposed", col="The Proposed"), size=0.75)+
                                        labs(y= ylab, x = xlab, sep="")+
                                        theme_bw()+

                                        theme(legend.position = c(.2,.75),text = element_text(size = 20),
                                              plot.title = element_blank()) +
                                        scale_linetype_manual(name = "",values = c("dotted","twodash"),
                                                              labels = c("Naive","The Proposed")) +
                                        scale_color_manual(name = "",values = c("#D95F02", "#1B9E77"),
                                                           labels = c("Naive","The Proposed")))

  }

  dat_graphs <- data.frame(years,
                           mean.My10_all,
                           se.My10_all,
                           CI_NE_U_all,CI_NE_L_all)


  #--Define title and axis labels:
  main <- paste("The CIF,","n=",sample_size, sep = "")
  xlab <- expression(paste(Age))
  ylab <- expression(paste(G(Age)))

  #--Line styles and colours:
  col <- c("#1B9E77", "#D95F02")
  lty <- 1:2	# Line type
  lwd <- 3.5	# Line width


  plot_list$Proposed <- print(ggplot(dat_graphs,aes(x = years)) +
                                geom_ribbon(aes(ymin=CI_NE_L_all,ymax=CI_NE_U_all),fill=col[1],alpha=0.3)+
                                #geom_ribbon(aes(ymin=CI_NA_L,ymax=CI_NA_U),fill=col[2],alpha=0.3)+

                                #geom_line(aes( y= mean.Nelson10,linetype ="Naive", col="Naive"), size=0.75)+
                                geom_line(aes( y= mean.My10_all,linetype = "The Proposed", col="The Proposed"), size=0.75)+
                                labs(y= ylab, x = xlab, sep="")+
                                theme_bw()+

                                theme(legend.position = c(.2,.75),text = element_text(size = 20),
                                      plot.title = element_blank()) +
                                scale_linetype_manual(name = "",values = c("twodash"),
                                                      labels = c("The Proposed")) +
                                scale_color_manual(name = "",values = c("#1B9E77"),
                                                   labels = c("The Proposed")))


  return(plot_list)
}



