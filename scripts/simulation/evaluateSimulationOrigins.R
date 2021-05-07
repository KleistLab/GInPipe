# Evaluation of Marias simulation of populations growth and estimating theta with the number of origins

# clear the environment variables
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(latex2exp)
library(ggplot2)
# load other r routines
source("splineRoutines.r")
source("readSimuData.r")

# to print float number as they are
options(scipen = 999)

#number of simulations
Nsim=10
#max size of bins
Nbins <- 12

#growthRate<-"0.52"
#growthRate_2 <-c("0.49","0.52")
growthRate<-"1.05"
growthRate_2 <-c("1.05","0.95")
switchLikeOrigin<- c("True", "False")

#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_allowMultipleMutatationPerSite"
#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_sinusRepl_2"
#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_testNewMLE"
#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_testOldMLEfirstOccurrence"
#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_testOldMLEfirstOccurrence_sinoid"
#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_testNewTotMLEallOccurrence_sinoid"

datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_testOldMLEallOccurrence_sinoid"
#datadir_all <- "/Users/msmith/Documents/RKI/nCov/data/simulation/simulation_effPopSize_testOldMLEallOccurrence"


#mutation_rate, and theta is given by 2*N*mu
mu=0.0001

# rescale for secons y axis
rescale <- function(x, minX, maxX, minY, maxY) {
  return((maxY-minY)*(x - minX)/(maxX-minX) + minY)
}


#Coefficient of quartile variation (cqv): (q3 - q1) / (q3 + q1)
get_cqv <- function(data) {
  data_qartiles <- as.vector(quantile(data, c(0.25, 0.5, 0.75)))
  return((data_qartiles[3] - data_qartiles[1]) / (data_qartiles[3] + data_qartiles[1]))
}

#Quartile based coefficient of variation (qcv): (q3 - q1)/q2
get_qcv <- function(data) {
  data_qartiles <- as.vector(quantile(data, c(0.25, 0.5, 0.75)))
  return((data_qartiles[3] - data_qartiles[1]) / data_qartiles[2])
}

branded_colors <- list(
  "blue"   = "#00798c",
  "red"    = "#d1495b",
  "yellow" = "#edae49",
  "green"  = "#66a182",
  "navy"   = "#2e4057", 
  "grey"   = "#8d96a3"
)

blues <- RColorBrewer::brewer.pal(4, "Blues")[c(2,4)]
greens <- RColorBrewer::brewer.pal(4, "Greens")[c(2,4)]
purples <-  RColorBrewer::brewer.pal(4, "Purples")[c(2,4)]

#for(sw in switchLikeOrigin) {
 sw='True'
  #for(gr2 in growthRate_2) {
   #gr2="0.95"
    gr2="1.05"
    #outpath_pub<-paste0(datadir_all, "/result_plots_publi_optiK_", gr2)
    outpath_pub<-paste0(datadir_all, "/result_plots_publi_optiK")
    datadir <- paste0(datadir_all,"/python_out_p_repl_",growthRate,"_repl2_",gr2,"_p_mut_",mu,"_switchInit_",sw)
    outputpath <- outpath_pub#paste0(datadir,"/result_plots_log_poi")
    
    dispersionFile <- paste0(outputpath, "/dispersion.csv")
    
    if(!dir.exists(outpath_pub))
      dir.create(outpath_pub, showWarnings = FALSE, recursive = TRUE)
    
    if(!dir.exists(outputpath))
      dir.create(outputpath, showWarnings = FALSE, recursive = TRUE)
    ##########################################################################
    ## Read Simulation Data for binning and subsampling 
    ##########################################################################
    simuTables = readData(datadir = datadir, Nsim = Nsim, Nbins= Nbins)
    
    all_simuData.table <- simuTables$allSimu
    all_splinePredictions.table <- simuTables$allSpline

    trueN_splinePrediction.table <- simuTables$trueSpline
    trueN_simuData.table <- simuTables$trueSimu
    
    minDate <- as.Date(min(trueN_simuData.table$meanBinDate))
    
    #read data for WS and WSMB with different subsamplings
    simuTables_sub <- readData_subsampling(datadir, Nsim, Nbins)
    all_simuData_sub.table <- simuTables_sub$allSimu
    all_splinePredictions_sub.table <- simuTables_sub$allSpline
    
    ##########################################################################
    ## Figures True Theta / N
    ##########################################################################
    
    #### TODO: Figure 1
    p_trueN <- ggplot() +
      geom_ribbon(data = trueN_splinePrediction.table, aes(x=t, ymin=lower, ymax=upper, fill=factor(simu)), col=NA, alpha=0.5, show.legend = F) +
      geom_line(data = trueN_splinePrediction.table, 
                #aes(x=days.as.Date(t,minDate), y=value, col=factor(simu)),
                aes(x=t, y=value, col=factor(simu)),
                size=1,
                show.legend = FALSE)+
      geom_point(aes(x=trueN_simuData.table$t,
                     y = trueN_simuData.table$trueN , 
                     col=factor(trueN_simuData.table$simu)), 
                 alpha=0.5,
                 show.legend = FALSE) +
      #scale_x_date(date_labels = "%B %Y")+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab("number of generations")+
      ylab(expression(N["true"]))
    
    #ggsave(filename = paste0(outputpath, "/trueN_splines_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_trueN, device = "pdf",width=24, height=15, units = "cm")
    ggsave(filename = paste0(outputpath, "/trueN_splines_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_trueN, device = "pdf",width=15, height=15, units = "cm")
    
    
    p_trueTheta <- ggplot() +
      geom_line(data = trueN_splinePrediction.table, 
                aes(x=t, y=value*2*mu, col=factor(simu)),
                size=1,
                show.legend = FALSE)+
      geom_point(data=trueN_simuData.table, aes(x=t,
                     y = trueN*2*mu, col=factor(simu)), alpha=0.5, show.legend = FALSE) +
      #scale_x_date(date_labels = "%B %Y")+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab("")+
      ylab(expression(theta["true"]))
  
    # ggsave(filename = paste0(outputpath, "/trueTheta_splines_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_trueTheta, device = "pdf",width=10, height=15, units = "cm")
    
    ##########################################################################
    ## Figures Theta estimation with all point estimates
    ##########################################################################

    NS_splinePredictions.table <- all_splinePredictions.table[all_splinePredictions.table$sample=="NS",]
    NS_simuData.table <- all_simuData.table[all_simuData.table$sample=="NS",]
    
    p_thetaNS <- ggplot(NS_splinePredictions.table, aes(x=t, y=value, linetype=sample, color=factor(simu))) +
      #geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=factor(simu)), col=NA, alpha=0.3, show.legend = F) +
      geom_line(size=1,show.legend = FALSE)+
      # geom_point(data=NS_simuData.table,aes(x=t, y=value),
      #            alpha=0.5,
      #            show.legend = FALSE)+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab("number of generations")+
      ylab(expression(theta["NS"]))
    
    ggsave(filename = paste0(outputpath, "/thetaNS_spline_noPoints_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_thetaNS, device = "pdf",width=10, height=15, units = "cm")
    
    ### true vs estimate
    p_Ntrue_vs_thetaNS <- ggplot(data=NS_simuData.table,aes(x=trueN_simuData.table$trueN, y=value))+
      geom_point(col="darkgreen", alpha=0.7, size=3)+
      geom_smooth(col="darkred", size=2, method='lm') +
      #ylim(0,7.5)+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab(expression(N["true"])) +
      ylab(expression(theta["NS"]))
    
    ggsave(filename = paste0(outputpath, "/thetaTrue_vs_thetaNS_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaNS, device = "pdf",width=10, height=10, units = "cm")
    
    p_Ntrue_vs_thetaNS_spline <- ggplot(data=NS_splinePredictions.table,aes(x=trueN_splinePrediction.table$value, y=value)) +
      geom_point(alpha=0.7, col="darkgreen", size=3)+
      geom_smooth(col="darkred", size=2, method='lm') +
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      ylim(layer_scales(p_Ntrue_vs_thetaNS)$y$range$range)+
      xlim(layer_scales(p_Ntrue_vs_thetaNS)$x$range$range)+
      xlab(expression(theta["true"])) +
      ylab(expression(theta["NS"]))
    
    ggsave(filename = paste0(outputpath, "/thetaTrue_vs_thetaNS_spline_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaNS_spline, device = "pdf",width=10, height=10, units = "cm")
  
    
    # Boxplot of quotient
    
   
    NS_simuData_quotient = NS_simuData.table$value/trueN_simuData.table$trueN
    NS_spline_quotient = NS_splinePredictions.table$value/trueN_splinePrediction.table$value
    
    #median quotient
    NS_median_simuData <- median(NS_simuData_quotient)
    NS_median_spline <- median(NS_spline_quotient)
    
    NS_simuData_cqv <- get_cqv(NS_simuData_quotient) 
    NS_simuData_qcv <- get_qcv(NS_simuData_quotient)

    NS_spline_cqv <- get_cqv(NS_spline_quotient) 
    NS_spline_qcv <- get_qcv(NS_spline_quotient)
    
    print(paste("NS QCV: ", format(NS_simuData_qcv, digits = 2), format(NS_spline_qcv, digits = 2)))
    print(paste("NS CQV: ", format(NS_simuData_cqv, digits = 2), format(NS_spline_cqv, digits = 2)))
    dispersionTable <- data.frame(sample="NS", QCV_sample = NS_simuData_qcv, QCV_fit = NS_spline_qcv,
                                  CQV_sample = NS_simuData_cqv, CQV_fit = NS_spline_cqv)
    
    NS_quotient.df <- data.frame(quotient=NS_simuData_quotient, data="point estimates")
    NS_quotient.df <- rbind(NS_quotient.df, data.frame(quotient=NS_spline_quotient, data="spline estimates"))
    
    p_log_quotient_bp_NS <-  ggplot(NS_quotient.df)+
      geom_violin(aes(x=data, y=log(quotient),fill=data), draw_quantiles = c(0.5), alpha = 0.5, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(NS_simuData_cqv, digits=2),'\n',"QCV: ", format(NS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(NS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(NS_spline_cqv, digits=2),'\n',"QCV: ", format(NS_spline_qcv, digits=2)),
      #              x=2.3, y=log(quantile(NS_spline_quotient, 0.9)) )) +
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      scale_fill_manual(values = greens)+
      xlab("") +
      ylab(expression(paste (log," (", theta["NS"]/N["true"], ")")))
    
    ggsave(filename = paste0(outputpath, "/log_NStheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_log_quotient_bp_NS, device = "pdf",width=15, height=15, units = "cm")
    
    p_quotient_bp_NS <-  ggplot(NS_quotient.df)+
      geom_violin(aes(x=data, y=(quotient),fill=data), draw_quantiles = c(0.5), alpha = 0.5, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(NS_simuData_cqv, digits=2),'\n',"QCV: ", format(NS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(NS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(NS_spline_cqv, digits=2),'\n',"QCV: ", format(NS_spline_qcv, digits=2)),
      #              x=2.3, y=log(quantile(NS_spline_quotient, 0.9)) )) +
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      scale_fill_manual(values = greens)+
      xlab("") +
      ylab(expression(paste (theta["NS"]/N["true"])))
    
    ggsave(filename = paste0(outputpath, "/NStheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_quotient_bp_NS, device = "pdf",width=15, height=15, units = "cm")
    
    
    # QQ Plot: Nötig? 
    # q_NS_1 <- quantile(NS_splinePredictions.table$value[NS_splinePredictions.table$simu==3], seq(0,1, 0.01))
    # q_true_1 <- quantile(trueN_splinePrediction.table$value[trueN_splinePrediction.table$simu==3], seq(0,1,0.01))
    # 
    # ggplot() +
    #   geom_point(aes(x=log(q_true_1), y=log(q_NS_1)))
    
    
    ##########################################################################
    ## Figures Theta estimation with subsampling
    ##########################################################################
  
    #plots woth random subsampling
    
    WS_splinePredictions.table <- all_splinePredictions.table[all_splinePredictions.table$sample=="WS",]
    WS_simuData.table <- all_simuData.table[all_simuData.table$sample=="WS",]
   
     p_thetaWS <- ggplot(WS_splinePredictions.table, 
                        aes(x=t, y=value, color=factor(simu))) +
      #geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=factor(simu)), col=NA, alpha=0.3, show.legend = F) +
      geom_line(size=1, show.legend = F)+
      # geom_point(data=WS_simuData.table,
      #            aes(x=t, y=value),
      #            alpha=0.5,
      #            show.legend = FALSE)+
      #scale_x_date(date_labels = "%B %Y")+
       theme_classic()+
       theme(axis.title = element_text(size = 20),
             axis.text.y = element_text(size=14),
             axis.text.x = element_text(size=14),
             legend.text = element_text(size=14),
             legend.title = element_text(size=14))+
      xlab("number of generations")+
      ylab(expression(theta["WS"]))
    
    
    ggsave(filename = paste0(outputpath, "/thetaWS_spline_noPoints_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_thetaWS, device = "pdf",width=10, height=15, units = "cm")
    
    
    # Boxplot of quotient
    
    WS_simuData_quotient = (WS_simuData.table$value)/(trueN_simuData.table$trueN)
    WS_spline_quotient = (WS_splinePredictions.table$value)/(trueN_splinePrediction.table$value)
    
    #median quotient
    WS_median_simuData <- median(WS_simuData_quotient)
    WS_median_spline <- median(WS_spline_quotient)
    
    
    WS_simuData_cqv <- get_cqv(WS_simuData_quotient) 
    WS_simuData_qcv <- get_qcv(WS_simuData_quotient)
    
    WS_spline_cqv <- get_cqv(WS_spline_quotient) 
    WS_spline_qcv <- get_qcv(WS_spline_quotient)
    
    print(paste("WS QCV: ", format(WS_simuData_qcv, digits = 2), format(WS_spline_qcv, digits = 2)))
    print(paste("WS CQV: ", format(WS_simuData_cqv, digits = 2), format(WS_spline_cqv, digits = 2)))
    dispersionTable <- rbind(dispersionTable, data.frame(sample="WS", QCV_sample = WS_simuData_qcv, QCV_fit = WS_spline_qcv,
                                  CQV_sample = WS_simuData_cqv, CQV_fit = WS_spline_cqv))
    
    WS_quotient.df <- data.frame(t=WS_simuData.table$t, quotient=WS_simuData_quotient, data="point estimates")
    WS_quotient.df <- rbind(WS_quotient.df, data.frame(t=WS_splinePredictions.table$t, quotient=WS_spline_quotient, data="spline estimates"))
    
    p_log_quotient_bp_WS <-  ggplot(WS_quotient.df[WS_quotient.df$quotient>=0,], aes(x=data, y=log(quotient)))+
      geom_violin(aes(fill=data), draw_quantiles = c(0.5), alpha = 0.7, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(WS_simuData_cqv, digits=2),'\n',"QCV: ", format(WS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(WS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(WS_spline_cqv, digits=2),'\n',"QCV: ", format(WS_spline_qcv, digits=2)),
      #               x=2.3, y=log(quantile(WS_spline_quotient, 0.9)) )) +
      theme_classic()+
      scale_colour_gradient(low = "yellow", high = "darkgreen")+ 
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      scale_fill_manual(values = blues)+
      xlab("") +
      ylab(expression(paste (log," (", theta["WS"]/N["true"], ")")))
    
    ggsave(filename = paste0(outputpath, "/log_WStheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_log_quotient_bp_WS, device = "pdf",width=15, height=20, units = "cm")
    
    p_quotient_bp_WS <-  ggplot(WS_quotient.df, aes(x=data, y=(quotient)))+
      geom_violin(aes(fill=data), draw_quantiles = c(0.5), alpha = 0.7, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(WS_simuData_cqv, digits=2),'\n',"QCV: ", format(WS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(WS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(WS_spline_cqv, digits=2),'\n',"QCV: ", format(WS_spline_qcv, digits=2)),
      #               x=2.3, y=log(quantile(WS_spline_quotient, 0.9)) )) +
      theme_classic()+
      scale_colour_gradient(low = "yellow", high = "darkgreen")+ 
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      scale_fill_manual(values = blues)+
      xlab("") +
      ylab(expression(paste (theta["WS"]/N["true"])))
    
    ggsave(filename = paste0(outputpath, "/WStheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_quotient_bp_WS, device = "pdf",width=15, height=20, units = "cm")
    
    
    # show jitter plot for the time
    p_quotient_jit_WS <- ggplot(WS_quotient.df[WS_quotient.df$quotient>0 ,], aes(x=data, y=log(quotient)))+
      geom_jitter(aes(col=t), width=0.2, alpha=0.2)+
      geom_violin(draw_quantiles = c(0.5), alpha = 0, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(WS_simuData_cqv, digits=2),'\n',"QCV: ", format(WS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(WS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(WS_spline_cqv, digits=2),'\n',"QCV: ", format(WS_spline_qcv, digits=2)),
      #               x=2.3, y=log(quantile(WS_spline_quotient, 0.9)) )) +
      theme_classic()+
      scale_colour_gradient2(low = "blue", mid= "yellow", high = "red", midpoint=max(WS_quotient.df$t)/2)+ 
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      #scale_fill_manual(values = blues)+
      xlab("") +
      ylab(expression(paste (log," (", theta["WS"]/N["true"], ")")))
    
    ggsave(filename = paste0(outputpath, "/log_WStheta_ratio_jitter_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_quotient_jit_WS, device = "pdf",width=15, height=15, units = "cm")
    
    # for the width of the point to be maximal 0.5
    maxNtrue_WS <- max(c(trueN_simuData.table$trueN, trueN_splinePrediction.table$value))
    width_WS <-((c(trueN_simuData.table$trueN, trueN_splinePrediction.table$value)[WS_quotient.df$quotient>0])/maxNtrue_WS)*0.48
    
    # show jitter plot for the size of the true value
    p_log_quotient_jit_val_WS <- ggplot(WS_quotient.df[WS_quotient.df$quotient>0,], aes(x=data, y=log(quotient)))+
      geom_jitter(aes(col=(c(trueN_simuData.table$trueN, trueN_splinePrediction.table$value)[WS_quotient.df$quotient>0])), width=width_WS, alpha=0.2)+
      geom_violin(draw_quantiles = c(0.5), alpha = 0, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(WS_simuData_cqv, digits=2),'\n',"QCV: ", format(WS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(WS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(WS_spline_cqv, digits=2),'\n',"QCV: ", format(WS_spline_qcv, digits=2)),
      #               x=2.3, y=log(quantile(WS_spline_quotient, 0.9)) )) +
      theme_classic()+
      scale_colour_gradient(low = "blue", high = "red")+ 
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      #scale_fill_manual(values = blues)+
      labs(col=expression(N["true"]))+
      xlab("") +
      ylab(expression(paste (log," (", theta["WS"]/N["true"], ")")))
    
    ggsave(filename = paste0(outputpath, "/log_WStheta_ratio_jitter_trueVal_repl2_",gr2,"_switchOrigin_",sw,"_width.pdf"), plot = p_log_quotient_jit_val_WS, device = "pdf",width=15, height=15, units = "cm")
    
    
    width_WS <-((c(trueN_simuData.table$trueN, trueN_splinePrediction.table$value))/maxNtrue_WS)*0.48
    
    p_quotient_jit_val_WS <- ggplot(WS_quotient.df, aes(x=data, y=(quotient)))+
      geom_jitter(aes(col=(c(trueN_simuData.table$trueN, trueN_splinePrediction.table$value))), width=width_WS, alpha=0.2)+
      geom_violin(draw_quantiles = c(0.5), alpha = 0, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(WS_simuData_cqv, digits=2),'\n',"QCV: ", format(WS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(WS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(WS_spline_cqv, digits=2),'\n',"QCV: ", format(WS_spline_qcv, digits=2)),
      #               x=2.3, y=log(quantile(WS_spline_quotient, 0.9)) )) +
      theme_classic()+
      scale_colour_gradient(low = "blue", high = "red")+ 
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      #scale_fill_manual(values = blues)+
      labs(col=expression(N["true"]))+
      xlab("") +
      ylab(expression(paste(theta["WS"]/N["true"])))
    
    ggsave(filename = paste0(outputpath, "/WStheta_ratio_jitter_trueVal_repl2_",gr2,"_switchOrigin_",sw,"_width.pdf"), plot = p_quotient_jit_val_WS, device = "pdf",width=15, height=15, units = "cm")
    
    
    
    p_quotient_jit_val_WS <- ggplot(WS_quotient.df[WS_quotient.df$quotient>0,], aes(x=data, y=log(quotient)))+
      geom_jitter(aes(col=(c(trueN_simuData.table$trueN, trueN_splinePrediction.table$value)[WS_quotient.df$quotient>0])), width=0.2, alpha=0.2)+
      geom_violin(draw_quantiles = c(0.5), alpha = 0, trim=T, show.legend = F)+
      #geom_label(aes(label = paste0("CQV: ", format(WS_simuData_cqv, digits=2),'\n',"QCV: ", format(WS_simuData_qcv, digits=2)),
      #               x=1.3, y=log(quantile(WS_simuData_quotient, 0.9)) )) +
      #geom_label(aes(label = paste0("CQV: ", format(WS_spline_cqv, digits=2),'\n',"QCV: ", format(WS_spline_qcv, digits=2)),
      #               x=2.3, y=log(quantile(WS_spline_quotient, 0.9)) )) +
      theme_classic()+
      scale_colour_gradient(low = "blue", high = "red")+ 
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      #scale_fill_manual(values = blues)+
      labs(col=expression(N["true"]))+
      xlab("") +
      ylab(expression(paste (log," (", theta["WS"]/N["true"], ")")))
    
    #ggsave(filename = paste0(outputpath, "/log_WStheta_ratio_jitter_trueVal_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_quotient_jit_val_WS, device = "pdf",width=15, height=15, units = "cm")
    
    
    
    ##### TODO! warum geht das nicht????
    #### subsampled data
    #ggplot() +
     # geom_line(data=all_splinePredictions.table[all_splinePredictions.table$sample=="WS",],# | all_splinePredictions.table$sample=="NS",], 
      #          aes(x=trueN_splinePrediction.table$value, y=value, col=factor(simu), linetype=factor(sample)), alpha=0.7, size=2)#, col="darkgreen")
    
    #### TODO: Figure 6 
    p_Ntrue_vs_thetaWS <-  ggplot()+
      geom_point(data=WS_simuData.table,aes(x=(trueN_simuData.table$trueN), y=(value), col=subsample),alpha=0.7, show.legend = F) +
      geom_smooth(data=WS_simuData.table,aes(x=trueN_simuData.table$trueN, y=value, weight=sampled_N), col="darkred", method=lm)+
     # geom_smooth(data=NS_simuData.table,aes(x=trueN_simuData.table$trueN, y=value), col="darkblue")+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            #legend.position="bottom",
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      scale_colour_gradient(low = "yellow", high = "darkgreen")+ 
      xlab(expression(N["true"])) +
      ylab(expression(theta["WS"]))
    
    ggsave(filename = paste0(outputpath, "/Ntrue_vs_thetaWS_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaWS, device = "pdf",width=10, height=12, units = "cm")
    
    notZero <- WS_simuData.table$value>0
    p_Ntrue_vs_thetaWS_log <-  ggplot(data=WS_simuData.table[notZero,],aes(x=log(trueN_simuData.table$trueN[notZero]), y=log(value), col=subsample))+
      geom_point(alpha=0.7) +
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      scale_colour_gradient(low = "yellow", high = "darkgreen")+ 
      xlab(expression(paste("log(", N["true"], ")"))) +
      ylab(expression(paste("log(", theta["WS"], ")")))
    #ggsave(filename = paste0(outputpath, "/Ntrue_vs_thetaWS_log_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaWS_log, device = "pdf",width=10, height=10, units = "cm")
    
    
    ####### plot results for different subsampling
    WS_simuData_sub.table <- all_simuData_sub.table[all_simuData_sub.table$sample=="WS",]
    trueTheta_WS <- trueN_simuData.table$trueN[unlist(sapply(seq(nrow(WS_simuData_sub.table)), function(i) which(trueN_simuData.table$t == WS_simuData_sub.table$t[i] & trueN_simuData.table$simu == WS_simuData_sub.table$simu[i])))]
    WS_splinePredictions_sub.table <- all_splinePredictions_sub.table[all_splinePredictions_sub.table$sample=="WS",]
    trueFit_WS <- trueN_splinePrediction.table$value[unlist(sapply(seq(nrow(WS_splinePredictions_sub.table)), function(i) which(trueN_splinePrediction.table$t == WS_splinePredictions_sub.table$t[i] & trueN_splinePrediction.table$simu == WS_splinePredictions_sub.table$simu[i])))]
    
    # Sieht nicht so gut aus, lieber für eine Beispiel zeigen wie sehr die einzelnen Abstufungen abweichen
    # ggplot()+
    #   # geom_point(data=WS_splinePredictions_sub.table,
    #   #            aes(x=rep(trueN_splinePrediction.table$value, length(unique(subsample))), y=value, col=factor(subsample)), alpha=0.7) +
    #   geom_point(data=WS_splinePredictions_sub.table[WS_splinePredictions_sub.table$subsample == 0.1,],
    #     aes(x=trueN_splinePrediction.table$value,  y=value,col=factor(subsample)), alpha=0.7)+
    #   geom_point(data=WS_splinePredictions_sub.table[WS_splinePredictions_sub.table$subsample == 0.5,],
    #              aes(x=rep(trueN_splinePrediction.table$value, length(unique(subsample))), y=value, col=factor(subsample)), alpha=0.7) +
    # geom_point(data=WS_splinePredictions_sub.table[WS_splinePredictions_sub.table$subsample == 0.9,],
    #            aes(x=rep(trueN_splinePrediction.table$value, length(unique(subsample))), y=value, col=factor(subsample)), alpha=0.7)
  
    
    simuWS=1
    # Plot as an exmpale for growth 1.05: 9,  growth=0.95: simulation 1, sinus: simulation 6/3
    p_WS_subsmpling <-  ggplot(WS_splinePredictions_sub.table[WS_splinePredictions_sub.table$simu==simuWS,]) +
      geom_line(data=NS_splinePredictions.table[NS_splinePredictions.table$simu==simuWS,], 
                aes(x=t, y=value, linetype=sample),size=2, show.legend = F)+
      geom_line(aes(x=t, y=value, alpha=subsample, group=factor(subsample)),
                size=2, linetype=2, color="darkred", show.legend = F)+
      #scale_x_date(date_labels = "%B %Y")+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            #legend.position="bottom",
            legend.title = element_text(size=14))+
      labs(linetype = element_blank())+
      xlab("")+
      ylab(expression(theta["WS"]))
    
    ggsave(filename = paste0(outputpath, "/thetaWS_subsampling_repl2_",gr2,"_switchOrigin_",sw,"_simu_",simuWS,".pdf"), plot = p_WS_subsmpling, device = "pdf",width=10, height=12, units = "cm")
    
    # Boxplot of quotients
    WS_simuData_quotient_sub = WS_simuData_sub.table$value/trueTheta_WS
    WS_spline_quotient_sub = WS_splinePredictions_sub.table$value/trueFit_WS
    
    WS_median_simuData_sub <- sapply(unique(WS_simuData_sub.table$subsample), function(s) median(WS_simuData_quotient_sub[which(WS_simuData_sub.table$subsample==s)]))
    WS_median_spline_sub <- sapply(unique(WS_splinePredictions_sub.table$subsample), function(s) median(WS_spline_quotient_sub[which(WS_splinePredictions_sub.table$subsample==s)]))
    
    WS_simuData_cqv_sub <- sapply(unique(WS_simuData_sub.table$subsample), function(s) get_cqv(WS_simuData_quotient_sub[which(WS_simuData_sub.table$subsample==s)])) 
    WS_simuData_qcv_sub <- sapply(unique(WS_simuData_sub.table$subsample), function(s) get_qcv(WS_simuData_quotient_sub[which(WS_simuData_sub.table$subsample==s)])) 
    
    WS_spline_cqv_sub <- sapply(unique(WS_splinePredictions_sub.table$subsample), function(s) get_cqv(WS_spline_quotient_sub[which(WS_splinePredictions_sub.table$subsample==s)])) 
    WS_spline_qcv_sub <- sapply(unique(WS_splinePredictions_sub.table$subsample), function(s) get_qcv(WS_spline_quotient_sub[which(WS_splinePredictions_sub.table$subsample==s)])) 
    
    print(paste("WS QCV: ", format(WS_simuData_qcv_sub, digits = 2), format(WS_spline_qcv_sub, digits = 2)))
    print(paste("WS CQV: ", format(WS_simuData_cqv_sub, digits = 2), format(WS_spline_cqv_sub, digits = 2)))
    
    WS_quotient_sub.df <- data.frame(t=WS_simuData_sub.table$t, subsample=WS_simuData_sub.table$subsample, quotient=WS_simuData_quotient_sub, data="point estimates")
    WS_quotient_sub.df <- rbind(WS_quotient_sub.df, data.frame(t=WS_splinePredictions_sub.table$t, subsample=WS_splinePredictions_sub.table$subsample, quotient=WS_spline_quotient_sub, data="spline estimates"))
    
    ggplot(WS_quotient_sub.df)+
      geom_boxplot(aes(x=subsample, y=log(quotient), fill=factor(data)))
    
    #violinplot
    p_quotient_bp_WS_sub <- ggplot(WS_quotient_sub.df, aes(x=subsample, y=log(quotient)))+
      geom_violin(aes(fill=factor(data)), draw_quantiles = c(0.5), alpha = 0.5, trim=T, position=position_dodge(width=0.7))+
      #geom_jitter(aes(col=t), width=0.25, alpha=0.1)
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_blank())+
      scale_fill_manual(values = blues)+
      xlab("Subsampling") +
      ylab(expression(paste (log," (", theta["WS"]/N["true"], ")")))
    
    ggsave(filename = paste0(outputpath, "/log_WStheta_sub_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_quotient_bp_WS_sub, device = "pdf",width=30, height=12, units = "cm")
   
    ##########################################################################
    ## Figures Theta estimation with binning
    ###########################################################################
    
    MB_splinePredictions.table <- all_splinePredictions.table[all_splinePredictions.table$sample=="MB",]
    MB_simuData.table <- all_simuData.table[all_simuData.table$sample=="MB",]
    
    p_thetaMB <- ggplot(MB_splinePredictions.table, 
                        aes(x=t, y=value, color=factor(simu))) +
      #geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=factor(simu)), col=NA, alpha=0.3, show.legend = F) +
      geom_line(size=1,show.legend = FALSE)+
      # geom_point(data=MB_simuData.table,
      #            aes(x=t, y=value, size=bin_size),
      #            alpha=0.5,
      #            show.legend = FALSE)+
      #scale_x_date(date_labels = "%B %Y")+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab("number of generations")+
      ylab(expression(theta["MB"]))
    
    ggsave(filename = paste0(outputpath, "/thetaMB_spline_noPoints_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_thetaMB, device = "pdf",width=10, height=15, units = "cm")
    
   
    #### TODO: Figure 4
    #true fitted values for the given bin time (as we do not have every intermediate data points)
    trueTheta_MB <- trueN_splinePrediction.table$value[unlist(sapply(seq(nrow(MB_simuData.table)), function(i) which(trueN_splinePrediction.table$t == MB_simuData.table$t[i] & trueN_splinePrediction.table$simu == MB_simuData.table$simu[i])))]
    trueFit_MB <- trueN_splinePrediction.table$value[unlist(sapply(seq(nrow(MB_splinePredictions.table)), function(i) which(trueN_splinePrediction.table$t == MB_splinePredictions.table$t[i] & trueN_splinePrediction.table$simu == MB_splinePredictions.table$simu[i])))]
    
    
    p_Ntrue_vs_thetaMB <- ggplot(data=MB_simuData.table, aes(x=trueTheta_MB, y=value))+
      geom_point(aes(size=factor(bin_size), col=factor(bin_size)), alpha=0.7) +
      geom_smooth(aes(weight=sampled_N),color="black", method="lm")+
      #geom_smooth(aes(color=factor(bin_size), fill=factor(bin_size)), linetype=1, alpha=0.1, show.legend = F)+
      #scale_colour_gradient(low = "white", high = "red")+ 
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab(expression(N["true"])) +
      ylab(expression(theta["MB"])) +
      labs(col="bin size", size="bin size")
    
    ggsave(filename = paste0(outputpath, "/thetaTrue_vs_thetaMB_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaMB, device = "pdf",width=10, height=10, units = "cm")
    
   # ratio plots
    MB_simuData_quotient = MB_simuData.table$value/trueTheta_MB
    MB_spline_quotient = MB_splinePredictions.table$value/trueFit_MB
    
    MB_median_simuData <- median(MB_simuData_quotient)
    MB_median_spline <- median(MB_spline_quotient)
    
    MB_simuData_cqv <- get_cqv(MB_simuData_quotient) 
    MB_simuData_qcv <- get_qcv(MB_simuData_quotient)
    
    MB_spline_cqv <- get_cqv(MB_spline_quotient) 
    MB_spline_qcv <- get_qcv(MB_spline_quotient)
    
    print(paste("MB QCV: ", format(MB_simuData_qcv, digits = 2), format(MB_spline_qcv, digits = 2)))
    print(paste("MB CQV: ", format(MB_simuData_cqv, digits = 2), format(MB_spline_cqv, digits = 2)))
    
    dispersionTable <- rbind(dispersionTable, data.frame(sample="MB", QCV_sample = MB_simuData_qcv, QCV_fit = MB_spline_qcv,
                                                         CQV_sample = MB_simuData_cqv, CQV_fit = MB_spline_cqv))
  
    
    
    p_log_boxplot_ratio_MB <- ggplot()+
      geom_violin(data=data.frame(data="point estimates", ratio=(MB_simuData_quotient), bin_size=factor(MB_simuData.table$bin_size)),
                  aes(x=data, y=log(ratio), fill=factor(bin_size)), draw_quantiles = c(0.5),alpha = 0.5)+
      geom_violin(data=data.frame(data="spline estimates", ratio=MB_spline_quotient), 
                  aes(x=data, y=log(ratio)), fill=blues[2], draw_quantiles = c(0.5), alpha = 0.5, show.legend = F)+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      ylab(expression(paste (log," (", theta["MB"]/N["true"], ")")))+
      xlab("")+
      labs(fill="bin size")
    ggsave(filename = paste0(outputpath, "/log_MBtheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_log_boxplot_ratio_MB, device = "pdf",width=20, height=10, units = "cm")
    
    p_boxplot_ratio_MB <- ggplot()+
      geom_violin(data=data.frame(data="point estimates", ratio=(MB_simuData_quotient), bin_size=factor(MB_simuData.table$bin_size)),
                  aes(x=data, y=(ratio), fill=factor(bin_size)), draw_quantiles = c(0.5),alpha = 0.5)+
      geom_violin(data=data.frame(data="spline estimates", ratio=MB_spline_quotient), 
                  aes(x=data, y=(ratio)), fill=blues[2], draw_quantiles = c(0.5), alpha = 0.5, show.legend = F)+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      ylab(expression(paste (theta["MB"]/N["true"])))+
      xlab("")+
      labs(fill="bin size")
    ggsave(filename = paste0(outputpath, "/MBtheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_boxplot_ratio_MB, device = "pdf",width=20, height=10, units = "cm")
    
       
    ##########################################################################
    ## Figures Theta estimation with subsampling and binning
    ##########################################################################
    WSMB_splinePredictions.table <- all_splinePredictions.table[all_splinePredictions.table$sample=="WSMB",]
    WSMB_simuData.table <- all_simuData.table[all_simuData.table$sample=="WSMB",]
    
    p_thetaWSMB <- ggplot(WSMB_splinePredictions.table, 
                        aes(x=t, y=value, linetype=sample, color=factor(simu))) +
      #geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=factor(simu)), col=NA, alpha=0.3, show.legend = F) +
      geom_line(size=1, show.legend = FALSE)+
      # geom_point(data=all_simuData.table[all_simuData.table$sample=="WSMB",],
      #            aes(x=t, y=value, size=bin_size),
      #            alpha=0.5,
      #            show.legend = FALSE)+
      #scale_x_date(date_labels = "%B %Y")+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab("number of generations")+
      ylab(expression(theta["WSMB"]))
    
    ggsave(filename = paste0(outputpath, "/thetaWSMB_spline_noPoints_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_thetaWSMB, device = "pdf",width=10, height=15, units = "cm")

       #### binning & subsampling
    
    ### for each estimate
    #true fitted values for the given bin time (as we do not have every intermediate data points)

    trueTheta_WSMB <- trueN_splinePrediction.table$value[unlist(sapply(seq(nrow(WSMB_simuData.table)), function(i) which(trueN_splinePrediction.table$t == WSMB_simuData.table$t[i] & trueN_splinePrediction.table$simu == WSMB_simuData.table$simu[i])))]
    trueFit_WSMB <- trueN_splinePrediction.table$value[unlist(sapply(seq(nrow(WSMB_splinePredictions.table)), function(i) which(trueN_splinePrediction.table$t == WSMB_splinePredictions.table$t[i] & trueN_splinePrediction.table$simu == WSMB_splinePredictions.table$simu[i])))]
    
    p_Ntrue_vs_thetaWSMB <- ggplot(data=WSMB_simuData.table, aes(x=trueTheta_WSMB, y=value))+
      geom_point(aes(size=factor(bin_size), col=factor(bin_size), alpha=subsample), show.legend = F) +
      geom_smooth(aes(weight=sampled_N),color="black", method="lm")+
      #geom_smooth(aes(color=factor(bin_size), fill=factor(bin_size)), linetype=1, alpha=0.1, show.legend = F)+
      #scale_colour_gradient(low = "white", high = "red")+ 
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      xlab(expression(N["true"])) +
      ylab(expression(theta["WSMB"])) +
      labs(col="bin size", size="bin size")
    
    ggsave(filename = paste0(outputpath, "/thetaTrue_vs_thetaWSMB_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaWSMB, device = "pdf",width=10, height=12, units = "cm")
    
    
      
    # ratio plots
    WSMB_simuData_quotient = WSMB_simuData.table$value/trueTheta_WSMB
    WSMB_spline_quotient = WSMB_splinePredictions.table$value/trueFit_WSMB
    
    WSMB_median_simuData <- median(WSMB_simuData_quotient)
    WSMB_median_spline <- median(WSMB_spline_quotient)
    
    WSMB_simuData_cqv <- get_cqv(WSMB_simuData_quotient) 
    WSMB_simuData_qcv <- get_qcv(WSMB_simuData_quotient)
    
    WSMB_spline_cqv <- get_cqv(WSMB_spline_quotient) 
    WSMB_spline_qcv <- get_qcv(WSMB_spline_quotient)
    
    print(paste("WSMB QCV: ", format(WSMB_simuData_qcv, digits = 2), format(WSMB_spline_qcv, digits = 2)))
    print(paste("WSMB CQV: ", format(WSMB_simuData_cqv, digits = 2), format(WSMB_spline_cqv, digits = 2)))
    dispersionTable <- rbind(dispersionTable, data.frame(sample="WSMB", QCV_sample = WSMB_simuData_qcv, QCV_fit = WSMB_spline_qcv,
                                                         CQV_sample = WSMB_simuData_cqv, CQV_fit = WSMB_spline_cqv))
    
    write.csv(dispersionTable, paste0(outputpath, "/dispersionTable_",gr2,".csv"), row.names = F)
    
    
    WSMB_ratio.df <- data.frame(data="point estimates", 
                          ratio=(WSMB_simuData_quotient), 
                          bin_size=factor(WSMB_simuData.table$bin_size),
                          value=WSMB_simuData.table$value, 
                          true=trueTheta_WSMB,
                          subsample=WSMB_simuData.table$subsample)
    
    WSMB_spline_ratio.df <- data.frame(data="spline estimates", ratio=WSMB_spline_quotient, true=trueFit_WSMB)
    
    # for the width of the point to be maximal 0.5
    maxNtrue_WSMB <- max(c(trueTheta_WSMB,trueFit_WSMB))
    width_point_WSMB <-(WSMB_ratio.df$true[WSMB_ratio.df$ratio != 0]/maxNtrue_WSMB)*0.5
    width_spline_WSMB <-(WSMB_spline_ratio.df$true[WSMB_spline_ratio.df$ratio != 0]/maxNtrue_WSMB)
    
    p_log_boxplot_ratio_WSMB <- ggplot()+
      geom_point(data=WSMB_ratio.df[WSMB_ratio.df$ratio != 0,],
                  aes(x=data, y=log(ratio), 
                      col=factor(bin_size), 
                      group=factor(bin_size), 
                      alpha=subsample,
                      size=true), 
                 position = position_jitterdodge(dodge.width = 0.9, jitter.width=width_point_WSMB))+
      geom_violin(data=WSMB_ratio.df[WSMB_ratio.df$ratio != 0,],
                  aes(x=data, y=log(ratio), group=factor(bin_size)), draw_quantiles = c(0.5),alpha = 0)+
      geom_point(data=WSMB_spline_ratio.df[WSMB_spline_ratio.df$ratio != 0,], 
                 aes(x=data, y=log(ratio),size=true/2), 
                 col=blues[2], alpha=0.5, shape=1,
                 position = position_jitterdodge(dodge.width = 0.9, jitter.width=width_spline_WSMB))+
      geom_violin(data=WSMB_spline_ratio.df[WSMB_spline_ratio.df$ratio != 0,], 
                  aes(x=data, y=log(ratio)), 
                  #fill=blues[2], 
                  draw_quantiles = c(0.5), alpha = 0, show.legend = F)+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      ylab(expression(paste (log," (", theta["WSMB"]/N["true"], ")")))+
      xlab("")+
      labs(col="bin size", size=expression(theta["true"]))
    ggsave(filename = paste0(outputpath, "/log_WSMBtheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_log_boxplot_ratio_WSMB, device = "pdf",width=30, height=15, units = "cm")
    
    
    width_point_WSMB <-(WSMB_ratio.df$true/maxNtrue_WSMB)*0.5
    width_spline_WSMB <-(WSMB_spline_ratio.df$true/maxNtrue_WSMB)
    
    p_boxplot_ratio_WSMB <- ggplot()+
      geom_point(data=WSMB_ratio.df,
                 aes(x=data, y=(ratio), 
                     col=factor(bin_size), 
                     group=factor(bin_size), 
                     alpha=subsample,
                     size=true), 
                 position = position_jitterdodge(dodge.width = 0.9, jitter.width=width_point_WSMB))+
      geom_violin(data=WSMB_ratio.df,
                  aes(x=data, y=(ratio), group=factor(bin_size)), draw_quantiles = c(0.5),alpha = 0)+
      geom_point(data=WSMB_spline_ratio.df[WSMB_spline_ratio.df$ratio != 0,], 
                 aes(x=data, y=(ratio),size=true/2), 
                 col=blues[2], alpha=0.5, shape=1,
                 position = position_jitterdodge(dodge.width = 0.9, jitter.width=width_spline_WSMB))+
      geom_violin(data=WSMB_spline_ratio.df[WSMB_spline_ratio.df$ratio != 0,], 
                  aes(x=data, y=(ratio)), 
                  #fill=blues[2], 
                  draw_quantiles = c(0.5), alpha = 0, show.legend = F)+
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14))+
      ylab(expression(paste (theta["WSMB"]/N["true"])))+
      xlab("")+
      labs(col="bin size", size=expression(theta["true"]))
    
    ggsave(filename = paste0(outputpath, "/WSMBtheta_ratio_boxplot_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_boxplot_ratio_WSMB, device = "pdf",width=30, height=15, units = "cm")
    
    
    
    
    
    # 
    # p_Ntrue_vs_thetaWSMB_log <-  ggplot(data=WSMB_simuData.table,
    #        aes(x=log(trueFit), y=log(value), col=factor(bin_size), size=factor(bin_size), alpha=subsample))+
    #   geom_point() +
    #   #scale_size(range = c(0, 10))+
    #   labs(col="bin size") +
    #   xlab(expression(paste("log(", theta["true"], ")"))) +
    #   ylab(expression(paste("log(", theta["WSMB"], ")")))
    # ggsave(filename = paste0(outputpath, "/trueTheta_vs_thetaWSMB_log_repl2_",gr2,"_switchOrigin_",sw,".pdf"), plot = p_Ntrue_vs_thetaWSMB_log, device = "pdf",width=24, height=15, units = "cm")
    # 

    # plot of different subsamplings
    data_WSMB_sub <- all_simuData_sub.table[all_simuData_sub.table$sample=="WSMB",]
    spline_WSMB_sub <- all_splinePredictions_sub.table[all_splinePredictions_sub.table$sample=="WSMB",]
    #true fitted values for the given bin time (as we do not have every intermediate data points)
    trueFit_sub <- trueN_splinePrediction.table$value[sapply(seq(nrow(data_WSMB_sub)), function(i) which(trueN_splinePrediction.table$t == data_WSMB_sub$t[i] & trueN_splinePrediction.table$simu == data_WSMB_sub$simu[i]))]
    trueFit_spline_sub <- trueN_splinePrediction.table$value[sapply(seq(nrow(spline_WSMB_sub)), function(i) which(trueN_splinePrediction.table$t == spline_WSMB_sub$t[i] & trueN_splinePrediction.table$simu == spline_WSMB_sub$simu[i]))]
    
    ggplot(spline_WSMB_sub) + 
      #geom_ribbon(aes(x=days.as.Date(t,minDate), ymin=value[subsample==0.1], ymax=value[subsample==0.9],color=factor(simu)))+
      geom_line(aes(x=days.as.Date(t,minDate), y=value, color=factor(simu), alpha=factor(subsample)), size=0.5)
    
   
    ##########################################################
    ###  Figures for presentation / publication 
    ##########################################################
    
    mycolors <- c("estid"="darkblue", "trueN"="darkred","esti"="blue", "true"="red")
  
    # scenario 1
    #s=2
    
    # scenario 2 1 +7
    #s=3
    
    # sinus: 1 
    s=1
    
    maxY1 <- max(NS_splinePredictions.table$value[NS_splinePredictions.table$simu==s])
    minY1 <- max(min(NS_splinePredictions.table$value[NS_splinePredictions.table$simu==s]),0)
    maxY2 <- max(trueN_simuData.table$trueN[trueN_simuData.table$simu==s])
    minY2 <- max(min(trueN_simuData.table$trueN[trueN_simuData.table$simu==s]),0)
    maxY2s <- max(trueN_splinePrediction.table$value[trueN_splinePrediction.table$simu==s])
    minY2s <- max(min(trueN_splinePrediction.table$value[trueN_splinePrediction.table$simu==s]),0)
    
    p_fig1_NS <- 
      ggplot() +
      # geom_ribbon(data=NS_splinePredictions.table[NS_splinePredictions.table$simu==s, ], 
      #               aes(x=days.as.Date(t,minDate), ymin =  sapply(lower, function(x) max(x, 0)), ymax = upper), fill = "grey70",alpha=0.5) +
      # # geom_line(data=trueN_splinePrediction.table[trueN_splinePrediction.table$simu==s,], aes(x=days.as.Date(t, minDate),
      #               y = rescale(value, minY2s, maxY2s, minY1, maxY1)))+
      geom_point(data=trueN_simuData.table[trueN_simuData.table$simu==s,], aes(x=t,
                     y = rescale(trueN, minY2s, maxY2s, minY1, maxY1)),
                     alpha=0.5,size=4, col=mycolors["trueN"],
                     show.legend = FALSE) +
      geom_line(data=NS_splinePredictions.table[NS_splinePredictions.table$simu==s, ], 
                aes(x=t, y=value),
                size=2, show.legend = F)+
      xlim(min(NS_splinePredictions.table$t), max(NS_splinePredictions.table$t))+
      #xlim(min(NS_splinePredictions.table$t), 96)+
      # scale_x_date(date_labels = "%B %Y", limits = c(days.as.Date(min(NS_splinePredictions.table$t), minDate), 
      #                                                days.as.Date(max(NS_splinePredictions.table$t), minDate)))+
      scale_y_continuous(
        expression(theta["est"]),
        sec.axis = sec_axis(~ rescale(., 0, maxY1, 0, maxY2s ), name = expression(N["true"]))) +
      theme_classic()+
      theme(axis.title = element_text(size = 20),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14),
            legend.text = element_text(size=14),
            legend.title = element_text(size=14),
            axis.title.y.right = element_text(color = mycolors["trueN"]),
            axis.text.y.right = element_text(color = mycolors["trueN"]))+
      xlab("number of generation")
    
    ggsave(filename = paste0(outpath_pub, "/fig1_tNS_v_nTrue_repl2_",gr2,"_sim_",s,".pdf"), plot = p_fig1_NS, device = "pdf",width=15, height=15, units = "cm")
    
    s=1
    
      maxY1_WSMB <- max(WSMB_splinePredictions.table$value[WSMB_splinePredictions.table$simu==s])
      minY1_WSMB <- max(min(WSMB_splinePredictions.table$value[WSMB_splinePredictions.table$simu==s]),0)
      maxY2 <- max(trueN_simuData.table$trueN[trueN_simuData.table$simu==s])
      minY2 <- max(min(trueN_simuData.table$trueN[trueN_simuData.table$simu==s]),0)
      maxY2s <- max(trueN_splinePrediction.table$value[trueN_splinePrediction.table$simu==s & trueN_splinePrediction.table$t<70])
      minY2s <- max(min(trueN_splinePrediction.table$value[trueN_splinePrediction.table$simu==s]),0)
      
      p_fig1_WSMB <- 
      ggplot() +
        # geom_ribbon(data=WSMB_splinePredictions.table[WSMB_splinePredictions.table$simu==s, ], 
        #             aes(x=days.as.Date(t,minDate), ymin = sapply(lower, function(x) max(x, 0)), ymax = upper), fill = "grey70",alpha=0.5) +
        # # geom_line(data=trueN_splinePrediction.table[trueN_splinePrediction.table$simu==s,], 
        #           aes(x=days.as.Date(t, minDate),
        #            y = rescale(value, minY2s, maxY2s, minY1_WSMB, maxY1_WSMB)), col="darkred")+
        geom_point(data=trueN_simuData.table[trueN_simuData.table$simu==s,], 
                   aes(x=t,
                       y = rescale(trueN, minY2s, maxY2s, minY1_WSMB, maxY1_WSMB)), 
                   alpha=0.5,size=4, col=mycolors["trueN"],
                   show.legend = FALSE) +
        geom_line(data=WSMB_splinePredictions.table[WSMB_splinePredictions.table$simu==s, ],
                  aes(x=t, y=value, linetype=sample),
                  size=2, show.legend = F)+
        xlim(min(WSMB_splinePredictions.table$t), max(WSMB_splinePredictions.table$t))+
        #xlim(min(WSMB_splinePredictions.table$t), 94)+
        # scale_x_date(date_labels = "%B %Y", limits = c(days.as.Date(min(WSMB_splinePredictions.table$t), minDate), 
        #                                                days.as.Date(max(WSMB_splinePredictions.table$t), minDate)))+
        scale_y_continuous(
          expression(theta["est"]),
          sec.axis = sec_axis(~ rescale(., 0, maxY1_WSMB, 0, maxY2s ), name = expression(N["true"]))) +
        theme_classic()+
        theme(axis.title = element_text(size = 20),
              axis.text.y = element_text(size=14),
              axis.text.x = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14),
              axis.title.y.right = element_text(color = mycolors["trueN"]),
              axis.text.y.right = element_text(color = mycolors["trueN"]))+
        xlab("number of generation")
      
      ggsave(filename = paste0(outpath_pub, "/fig1_tWSMB_nTrue_repl2_",gr2,"_sim_",s,".pdf"), plot = p_fig1_WSMB, device = "pdf",width=15, height=15, units = "cm")
      
    
      ### Figure 2: THere is a linear relationship
      
      p_fig2a_Ntrue_vs_thetaNS <- ggplot(data=NS_simuData.table,aes(x=trueN_simuData.table$trueN, y=value))+
        #p_fig2a_Ntrue_vs_thetaWSMB <- ggplot(data=WSMB_simuData.table,aes(x=trueTheta_WSMB, y=value))+
        geom_point(col="darkgreen", alpha=0.7, size=3)+
        geom_smooth(col="darkred", size=2) +
        #ylim(0,7.5)+
        theme_classic()+
        theme(axis.title = element_text(size = 20),
              axis.text.y = element_text(size=14),
              axis.text.x = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14))+
        xlab(expression(N["true"])) +
        ylab(expression(theta["est"]))
      ggsave(filename = paste0(outpath_pub, "/fig2a_NTrue_vs_tNS_repl2_",gr2,".pdf"), plot = p_fig2a_Ntrue_vs_thetaNS, device = "pdf",width=10, height=10, units = "cm")
      
      
      p_fig2b_Ntrue_thetaNS_log_ratio_violin<- ggplot(NS_quotient.df)+
        geom_violin(aes(x=data, y=log(quotient),fill=data), draw_quantiles = c(0.5), alpha = 0.5, trim=T, show.legend = F)+
        #geom_label(aes(label = paste0("CQV: ", format(NS_simuData_cqv, digits=2),'\n',"QCV: ", format(NS_simuData_qcv, digits=2)),
        #               x=1.3, y=log(quantile(NS_simuData_quotient, 0.9)) )) +
        #geom_label(aes(label = paste0("CQV: ", format(NS_spline_cqv, digits=2),'\n',"QCV: ", format(NS_spline_qcv, digits=2)),
        #              x=2.3, y=log(quantile(NS_spline_quotient, 0.9)) )) +
        theme_classic()+
        theme(axis.title = element_text(size = 20),
              axis.text.y = element_text(size=14),
              axis.text.x = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14))+
        scale_fill_manual(values = greens)+
        xlab("") +
        ylab(expression(paste (log," (", theta["est"]/N["true"], ")")))
      ggsave(filename = paste0(outpath_pub, "/fig2b_log_NSratio_violin_repl2_",gr2,".pdf"), plot = p_fig2b_Ntrue_thetaNS_log_ratio_violin, device = "pdf",width=15, height=15, units = "cm")
      
      p_fig2b_Ntrue_thetaNS_ratio_violin<- ggplot(NS_quotient.df)+
        geom_boxplot(aes(x=data, y=(quotient),fill=data), draw_quantiles = c(0.5), alpha = 0.5, trim=T, show.legend = F)+
        #geom_label(aes(label = paste0("CQV: ", format(NS_simuData_cqv, digits=2),'\n',"QCV: ", format(NS_simuData_qcv, digits=2)),
        #               x=1.3, y=log(quantile(NS_simuData_quotient, 0.9)) )) +
        #geom_label(aes(label = paste0("CQV: ", format(NS_spline_cqv, digits=2),'\n',"QCV: ", format(NS_spline_qcv, digits=2)),
        #              x=2.3, y=log(quantile(NS_spline_quotient, 0.9)) )) +
        theme_classic()+
        theme(axis.title = element_text(size = 20),
              axis.text.y = element_text(size=14),
              axis.text.x = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14))+
        scale_fill_manual(values = greens)+
        xlab("") +
        ylab(expression(paste (theta["est"]/N["true"])))
      ggsave(filename = paste0(outpath_pub, "/fig2b_NSratio_boxplot_repl2_",gr2,".pdf"), plot = p_fig2b_Ntrue_thetaNS_ratio_violin, device = "pdf",width=15, height=15, units = "cm")
      
      
      # NS_quotient.df_sc2 <- NS_quotient.df
      # WS_quotient.df_sc2 <- WS_quotient.df
      # WSMB_spline_ratio.df_sc2 <- WSMB_spline_ratio.df
      # trueN_splinePrediction.table_sc2 <- trueN_splinePrediction.table
      # trueFit_WSMB_sc2 <- trueFit_WSMB
      
      # NS_quotient.df_sc1 <- NS_quotient.df
      # WS_quotient.df_sc1 <- WS_quotient.df
      # WSMB_spline_ratio.df_sc1 <- WSMB_spline_ratio.df
      # trueFit_WSMB_sc1 <- trueFit_WSMB
      
      # for the width of the point to be maximal 0.5
      maxNtrue<- max(trueN_splinePrediction.table$value)

      ratio.df <- data.frame(quotient=WS_quotient.df$quotient[WS_quotient.df$data=="spline estimates"], trueN=trueN_splinePrediction.table$value, fill="scenario2", sample="subsampling")
     # ratio.df <- rbind(ratio.df, data.frame(quotient=NS_quotient.df$quotient[NS_quotient.df$data=="spline estimates"], trueN=trueN_splinePrediction.table$value, fill="scenario2", sample="all"))
      #ratio.df <- rbind(ratio.df, data.frame(quotient=NS_quotient.df_sc1$quotient[NS_quotient.df_sc1$data=="spline estimates"], fill="scenario1", sample="all"))
      #ratio.df <- rbind(ratio.df, data.frame(quotient=WS_quotient.df_sc1$quotient[WS_quotient.df_sc1$data=="spline estimates"], fill="scenario1", sample="subsampling"))
      ratio.df <- rbind(ratio.df, data.frame(quotient=WSMB_spline_ratio.df$ratio[WSMB_spline_ratio.df$data=="spline estimates"], trueN=trueFit_WSMB,fill="scenario2", sample="subsampling + binning"))
      #ratio.df <- rbind(ratio.df, data.frame(quotient=WSMB_spline_ratio.df_sc1$ratio[WSMB_spline_ratio.df_sc1$data=="spline estimates"], fill="scenario1", sample="subsampling + binning"))
      
      width_all <- ratio.df$trueN/maxNtrue*0.4
      
      p_fig3_Ntrue_thetas_ratio_violin <- 
        ggplot(data=ratio.df[ratio.df$quotient>=0,]) +
        geom_violin(aes(x=sample, y=(quotient),fill=fill), 
          position="dodge",draw_quantiles = c(0.5), alpha = 0.5, trim=T, show.legend = F)+
        theme_classic()+
        theme(axis.title = element_text(size = 20),
              axis.text.y = element_text(size=14),
              axis.text.x = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14))+
        scale_fill_manual(values = blues[2])+
        labs(x="", y=expression(paste (theta["est"]/N["true"])), fill="")

      ggsave(filename = paste0(outpath_pub, "/fig3_ratios_violin_WS",".pdf"), plot = p_fig3_Ntrue_thetas_ratio_violin, device = "pdf",width=30, height=15, units = "cm")
      
      
      p_fig3_Ntrue_thetas_ratio_violin_dots <- 
        ggplot(data=ratio.df[ratio.df$quotient>=0,]) +
        geom_jitter(aes(x=sample, y=(quotient),col=trueN, size=trueN), alpha=0.2, width=width_all)+
        geom_violin(aes(x=sample, y=(quotient),fill=fill), 
                    position="dodge",draw_quantiles = c(0.5), alpha = 0, trim=T, show.legend = F)+
        theme_classic()+
        scale_colour_gradient(low = "blue", high = "red")+ 
        #scale_colour_gradient2(low = "blue", mid="yellow", high = "red",  midpoint=max(trueFit_WSMB)/2)+ 
        theme(axis.title = element_text(size = 20),
              axis.text.y = element_text(size=14),
              axis.text.x = element_text(size=14),
              legend.text = element_text(size=14),
              legend.title = element_text(size=14))+
        scale_fill_manual(values = blues[2])+
        labs(x="", y=expression(paste (theta["est"]/N["true"])), fill="")
      
      ggsave(filename = paste0(outpath_pub, "/fig3_ratios_violin_dots_WS",".pdf"), plot = p_fig3_Ntrue_thetas_ratio_violin_dots, device = "pdf",width=25, height=15, units = "cm")
      
      
    
  #}
 #}     
        
      
      
   