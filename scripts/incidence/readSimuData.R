library(latex2exp)
library(ggplot2)
library(reshape2)
# load other r routines
source("splineRoutines.r")


readData <- function(datadir, Nsim, Nbins) {
  #collecting all data
  all_simuData.table <- data.frame()
  all_splinePredictions.table <- data.frame()
  trueN_splinePrediction.table <- data.frame()
  trueN_simuData.table <- data.frame()
  
  for(i in 1:Nsim) {
    #table with taking the complete sample
    input.table_NS = data.frame()
    file=paste0(datadir,"/NS_theta_origins_sim_",i,".tsv")
    if(file.exists(file)) {
      input.table_NS <- read.table(file,
                                   header = T, sep="\t")
      input.table_NS =  input.table_NS[-seq(4), ]
    }
    
    #table with taking a subsample
    input.table_WS = data.frame()
    file=paste0(datadir,"/WS_theta_origins_sim_",i,".tsv")
    if(file.exists(file)) {
      input.table_WS <- read.table(file,
                                   header = T, sep="\t")
      input.table_WS<- input.table_WS[-seq(4),]
    }
    
    #table with binning
    input.table_MB = data.frame()
    #table with binning and subsampling
    input.table_WSMB = data.frame()
    for(nb in seq(2,Nbins,2)) {
      file=paste0(datadir,"/MB_theta_origins_sim_",i,"_numBin_",nb,".tsv")
      if(file.exists(file)) {
        input.table_MB_tmp <- read.table(file,
                                         header = T, sep="\t")
        input.table_MB <- rbind(input.table_MB, input.table_MB_tmp[-seq(4),])
      }
      
      
      file=paste0(datadir,"/WSMB_theta_origins_sim_",i,"_numBin_",nb,".tsv")
      if(file.exists(file)) {
        input.table_WSMB_tmp <- read.table(file,
                                           header = T, sep="\t")
        input.table_WSMB<- rbind(input.table_WSMB, input.table_WSMB_tmp[-seq(4),])
      }
    }
    
    
    
    all_simuData.table = rbind(all_simuData.table, data.frame(data.frame(input.table_NS, sample="NS", simu=i, subsample=1)), 
                               data.frame(input.table_WS, sample="WS",simu=i, subsample=input.table_WS$sampled_N/input.table_WS$trueN), 
                               data.frame(input.table_MB, sample="MB",simu=i, subsample=1),
                               data.frame(input.table_WSMB, sample="WSMB",simu=i, subsample=input.table_WSMB$sampled_N/input.table_WSMB$trueN))
    #instead of true N compare true theta
    #trueN_simuData.table = rbind(trueN_simuData.table, data.frame(t=input.table_NS$t, true_N=input.table_NS$trueN*2*mu, simu=i))
    trueN_simuData.table = rbind(trueN_simuData.table, data.frame(t=input.table_NS$t, trueN=input.table_NS$trueN*2*mu, meanBinDate=input.table_NS$meanBinDate, simu=i))
    
    # Spline computations
    
    # vector for which the outcomes should be predicted
    #xx<- seq(min(input.table_NS$t), max(input.table_NS$t), len = max(input.table_NS$t) - min(input.table_NS$t)+1)
    
    #Splines for true N
    #trueN_spline <- computeSpline(data.frame(t=input.table_NS$t, value=input.table_NS$true_N*2*mu, variance=1))
    # trueN_spline <- computeSplineTable(data.frame(meanBinDate=input.table_NS$meanBinDate, value=input.table_NS$trueN, variance=1))
    # splinePred_trueN <- predict.gam(trueN_spline, data.frame(t=xx))
    # splinePred_deriv_trueN <- computeSplineDerivativeTable(xx, trueN_spline)
    splinePred_trueN <- computeSplineTable(data.frame(t=input.table_NS$t, meanBinDate=input.table_NS$meanBinDate, value=input.table_NS$trueN*2*mu, variance=1))
    
    #spline for all theta with all samples
    # simu_NS_spline <- computeSpline(input.table_NS)
    # splinePred_simu_NS <- predict.gam(simu_NS_spline, data.frame(t=xx))
    # splinePred_deriv_simu_NS <- computeSplineDerivativeTable(xx, simu_NS_spline)
    splinePred_simu_NS <-  computeSplineTable(input.table_NS)
    
    
    #spline for all theta with subsamples
    # simu_WS_spline <- computeSpline(input.table_WS)
    # splinePred_simu_WS <- predict.gam(simu_WS_spline, data.frame(t=xx))
    # splinePred_deriv_simu_WS <- computeSplineDerivativeTable(xx, simu_WS_spline)
    splinePred_simu_WS <- computeSplineTable(input.table_WS)
    
    #spline for all theta with binning
    # xx_MB <- xx[which(xx>min(input.table_MB$t))]
    # simu_MB_spline <- computeSpline(input.table_MB)
    # splinePred_simu_MB <- predict.gam(simu_MB_spline, data.frame(t=xx_MB))
    # splinePred_deriv_simu_MB <- computeSplineDerivativeTable(xx_MB, simu_MB_spline)
    splinePred_simu_MB <-  computeSplineTable(input.table_MB)
    
    #spline for all theta with subsamples and binning
    # xx_WSMB <- xx[which(xx>min(input.table_WSMB$t))]
    # simu_WSMB_spline <- computeSpline(input.table_WSMB)
    # splinePred_simu_WSMB <- predict.gam(simu_WSMB_spline, data.frame(t=xx_WSMB))
    # splinePred_deriv_simu_WSMB <- computeSplineDerivativeTable(xx_WSMB, simu_WSMB_spline)
    splinePred_simu_WSMB <- computeSplineTable(input.table_WSMB)
    
    # trueN_splinePrediction.table <- rbind(trueN_splinePrediction.table,
    #                                       data.frame(t=xx, splinePred=splinePred_trueN,splineDeriv= splinePred_deriv_trueN, sample="true N", simu=i)
    #                                       )
    trueN_splinePrediction.table <- rbind(trueN_splinePrediction.table,
                                          data.frame(splinePred_trueN, sample="true N", simu=i))
    
    all_splinePredictions.table <- rbind(all_splinePredictions.table, 
                                         data.frame(splinePred_simu_NS, sample="NS", simu=i), 
                                         data.frame(splinePred_simu_WS, sample="WS", simu=i), 
                                         data.frame(splinePred_simu_MB, sample="MB", simu=i),
                                         data.frame(splinePred_simu_WSMB, sample="WSMB", simu=i) 
    ) 
  }
 return(list(allSimu=all_simuData.table, trueSimu=trueN_simuData.table, allSpline=all_splinePredictions.table, trueSpline=trueN_splinePrediction.table)) 
}

readData_intros <- function(datadir, Nsim, Nbins, num_intros) {
  #collecting all data
  all_simuData.table <- data.frame()
  all_splinePredictions.table <- data.frame()
  trueN_splinePrediction.table <- data.frame()
  trueN_simuData.table <- data.frame()
  
  
  for(i in 1:Nsim) {
    
    # input.table_NS = data.frame()
    # file=paste0(datadir,"/NS_theta_origins_sim_",i,".tsv")
    # if(file.exists(file)) {
    #   input.table_NS <- read.table(file,header = T, sep="\t")
    #   input.table_NS =  input.table_NS[-seq(4), ]
    # }
    # 
    # input.table_WSMB = data.frame()
    # for(nb in seq(2,Nbins,2)) {
    #   file=paste0(datadir,"/WSMB_theta_origins_sim_",i,"_numBin_",nb,".tsv")
    #   input.table_WSMB_tmp <- read.table(file, header = T, sep="\t")
    #   input.table_WSMB<- rbind(input.table_WSMB, input.table_WSMB_tmp[-seq(4),])
    # }
    # all_simuData.table = rbind(all_simuData.table,
    #                            data.frame(input.table_WSMB, numIntros=0, totalOrigins=sum(input.table_WSMB$origins),
    #                                       sample="WSMB",simu=i, subsample=input.table_WSMB$sampled_N/input.table_WSMB$trueN, intro=0),
    #                            data.frame(input.table_NS, numIntros=0, totalOrigins=sum(input.table_NS$origins),
    #                                      sample="NS",simu=i, subsample=input.table_NS$sampled_N/input.table_NS$trueN, intro=0))
    # 
    # splinePred_simu_WSMB <- computeSplineTable(input.table_WSMB)
    # splinePred_simu_NS <- computeSplineTable(input.table_NS)
    # all_splinePredictions.table <- rbind(all_splinePredictions.table,
    #                                      data.frame(splinePred_simu_WSMB, sample="WSMB", simu=i, intro=0),
    #                                      data.frame(splinePred_simu_NS, sample="NS", simu=i, intro=0))
    # 
    # trueN_simuData.table = rbind(trueN_simuData.table, 
    #                              data.frame(t=input.table_NS$t, trueN=input.table_NS$trueN, meanBinDate=input.table_NS$meanBinDate, simu=i, intro=0))
    
    for(ni in num_intros) {
      input.table_NS_intro = data.frame()
      file=paste0(datadir,"/NS_theta_origins_sim_",i,"_intro_",ni,".tsv")
      if(file.exists(file)) {
        input.table_NS_intro <- read.table(file,
                                     header = T, sep="\t")
        input.table_NS_intro =  input.table_NS_intro[-seq(4), ]
      }
      
      #table with binning and subsampling
      # input.table_WSMB = data.frame()
      input.table_WSMB_intro = data.frame()
      for(nb in seq(2,Nbins,2)) {

        # input.table_WSMB_tmp <- read.table(paste0(datadir,"/WSMB_theta_origins_sim_",i,"_numBin_",nb,".tsv"),
        #                                    header = T, sep="\t")
        # input.table_WSMB<- rbind(input.table_WSMB, input.table_WSMB_tmp[-seq(4),])
        
        #with intros
        file=paste0(datadir,"/WSMB_theta_origins_sim_",i,"_numBin_",nb,"_intro_",ni,".tsv")
        if(file.exists(file)) {
          input.table_WSMB_tmp <- read.table(file,
                                             header = T, sep="\t")
          input.table_WSMB_intro<- rbind(input.table_WSMB_intro, input.table_WSMB_tmp[-seq(4),])
        }
      }
      
      all_simuData.table = rbind(all_simuData.table, 
                                 data.frame(input.table_NS_intro, totalOrigins=sum(input.table_NS_intro$origins), sample="NS",simu=i, subsample=1, intro=ni), 
                                 data.frame(input.table_WSMB_intro, totalOrigins=sum(input.table_WSMB_intro$origins), sample="WSMB",simu=i, subsample=input.table_WSMB_intro$sampled_N/input.table_WSMB_intro$trueN, intro=ni))
      #instead of true N compare true theta
      trueN_simuData.table = rbind(trueN_simuData.table, 
                                   #data.frame(t=input.table_NS$t, trueN=input.table_NS$trueN, meanBinDate=input.table_NS$meanBinDate, simu=i, intro=ni),
                                   data.frame(t=input.table_NS_intro$t, trueN=input.table_NS_intro$trueN*2*mu, meanBinDate=input.table_NS_intro$meanBinDate, simu=i, intro=ni))
      
      #splinePred_trueN <- computeSplineTable(data.frame(t=input.table_NS$t, meanBinDate=input.table_NS$meanBinDate, value=input.table_NS$trueN, variance=1))
      splinePred_trueN_intro <- computeSplineTable(data.frame(t=input.table_NS_intro$t, meanBinDate=input.table_NS_intro$meanBinDate, value=input.table_NS_intro$trueN*2*mu, variance=1))
      #splinePred_simu_NS <-  computeSplineTable(input.table_NS)
      splinePred_simu_NS_intro <-  computeSplineTable(input.table_NS_intro)
      #splinePred_simu_WSMB <- computeSplineTable(input.table_WSMB)
      splinePred_simu_WSMB_intro <- computeSplineTable(input.table_WSMB_intro)
      
      trueN_splinePrediction.table <- rbind(trueN_splinePrediction.table,
                                            #data.frame(splinePred_trueN, sample="true N", simu=i, intro=ni),
                                            data.frame(splinePred_trueN_intro, sample="true N", simu=i, intro=ni))
      
      all_splinePredictions.table <- rbind(all_splinePredictions.table, 
                                           #data.frame(splinePred_simu_NS, sample="NS", simu=i, intro=ni), 
                                           data.frame(splinePred_simu_NS_intro, sample="NS", simu=i, intro=ni), 
                                           #data.frame(splinePred_simu_WSMB, sample="WSMB", simu=i, intro=ni),
                                           data.frame(splinePred_simu_WSMB_intro, sample="WSMB", simu=i, intro=ni))
    }
  }
  return(list(allSimu=all_simuData.table, trueSimu=trueN_simuData.table, allSpline=all_splinePredictions.table, trueSpline=trueN_splinePrediction.table)) 
}

readData_subsampling <- function(datadir, Nsim, Nbins) {
  #collecting all data
  all_simuData.table <- data.frame()
  all_splinePredictions.table <- data.frame()
  
  for(i in 1:Nsim) {
    for(ss in as.character(seq(0.1, 0.9, 0.1))) {
     #ss=0.7 
      input.table_WS_sub = data.frame()
      file=paste0(datadir,"/WS_theta_origins_sim_",i,"_subsample_",ss,".tsv")
      if(file.exists(file)) {
        input.table_WS_sub <- read.table(file,header = T, sep="\t")
        #skip first 4 values, as the beginning is often too high
        input.table_WS_sub =  input.table_WS_sub[-seq(4), ]
      }
      
      #table with binning and subsampling
      input.table_WSMB_sub = data.frame()
      for(nb in seq(2,Nbins,2)) {
        file=paste0(datadir,"/WSMB_theta_origins_sim_",i,"_numBin_",nb,"_subsample_",ss,".tsv")
        if(file.exists(file)) {
          input.table_WSMB_tmp <- read.table(file,header = T, sep="\t")
          input.table_WSMB_sub<- rbind(input.table_WSMB_sub, input.table_WSMB_tmp[-seq(4),])
        }
      }
      
      all_simuData.table = rbind(all_simuData.table, 
                                 data.frame(input.table_WS_sub, totalOrigins=sum(input.table_WS_sub$origins), sample="WS",simu=i, subsample=ss), 
                                 data.frame(input.table_WSMB_sub, totalOrigins=sum(input.table_WSMB_sub$origins), sample="WSMB",simu=i, subsample=ss))
      
      splinePred_simu_WS_sub <-  computeSplineTable(input.table_WS_sub)
      splinePred_simu_WSMB_sub <- computeSplineTable(input.table_WSMB_sub)
      
      all_splinePredictions.table <- rbind(all_splinePredictions.table, 
                                           data.frame(splinePred_simu_WS_sub, sample="WS", simu=i, subsample=ss), 
                                           data.frame(splinePred_simu_WSMB_sub, sample="WSMB", simu=i, subsample=ss))
    }
  }
  return(list(allSimu=all_simuData.table, allSpline=all_splinePredictions.table)) 
}

