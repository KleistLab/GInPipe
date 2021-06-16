library(ggplot2)
library(scales)
#' Plot Re estimates.
#'
#' Plotting the Re estimtes with its confidence interval and saving it at the given location.
#'
#' @param re.table a data frame containing a column "value" with R_e estimates, values "lower" and "upper" for the confidence interval and a column "date" 
#' @param outputFile (optional) The Re plot is written to this file
#' @param title (optional) a vector of time points which should be interpolated


plot_Re <- function(re.table, outputFile=NULL, title=NULL){
  p_interp_d1 <- ggplot() +
    geom_line(aes(x=as.Date(re.table$date,"%Y-%m-%d"), y=re.table$value), colour="darkred", size=2, alpha=0.7)+
    geom_line(aes(x=as.Date(re.table$date,"%Y-%m-%d"), y=re.table$lower), size=0.5, alpha=0.7, linetype="dashed")+
    geom_line(aes(x=as.Date(re.table$date,"%Y-%m-%d"), y=re.table$upper), size=0.5, alpha=0.7, linetype="dashed")+
    #geom_point(aes(tt.df[,1], d1_gam), alpha=0.5, size=2) +
    #ylim(c(0,100))+
    xlab("") +
    scale_x_date(date_breaks = "2 months" , date_labels = "%b %Y") +
    ylab(expression(paste(R[0]))) +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  if(!is.null(title))
    p_interp_d1 <- p_interp_d1 + ggtitle(title)

  if(!is.null(outputFile))
    ggsave(p_interp_d1, height = 6, width = 10,dpi = 220, file = outputFile)
}


#' Plot smoothed phi estimates over time.
#'
#' Plotting the smoothed phi estimates, both with and without point estimates (if the output files are given). 
#' If present, the reported cases are added to plot with the scale on the right y-axis.
#'
#' @param smooth.table a data frame containing the smoothed Medians as output from smooth_point_estimates
#' @param point.table (optional) a data frame containing the point estimates for phi in the column "value", a "date", and the "sampleSize"
#' @param outputFile (optional) The plot without point estimates is written to this file
#' @param outputFileDots (optional) The plot with point estimates is written to this file
#' @param cases.table (optional) The data frame containing the reported cases with the rolling average in column "new_cases_avrg" and "date
#' @param title (optional) For the title of the plot.

plot_smoothed_phi_with_new_cases <- function(smooth.table, point.table=NULL, cases.table=NULL, outputFile=NULL, outputFileDots=NULL, title=NULL) {
    
  mycolors <- c("estid"="dodgerblue4", "trued"="darkred","esti"="dodgerblue4", "true"="red")
  
  p_smooth_esti_realN <- ggplot() +
    geom_line(aes(x=as.Date(smooth.table$date,"%Y-%m-%d"), y=smooth.table$smoothMedian), size=2.3, colour=mycolors["esti"], alpha=0.5)+
    geom_line(aes(x=as.Date(smooth.table$date,"%Y-%m-%d"), y=smooth.table$smooth5), size=1.3, colour=mycolors["esti"], alpha=0.55, linetype="dashed")+
    geom_line(aes(x=as.Date(smooth.table$date,"%Y-%m-%d"), y=smooth.table$smooth95), size=1.3, colour=mycolors["esti"], alpha=0.55, linetype="dashed")+
    xlab("")+
    scale_x_date(date_breaks = "1 months", date_labels = "%b %Y") +
    theme(
      axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
      axis.title.y = element_text(color = mycolors["estid"]),
      axis.text.y = element_text(color = mycolors["estid"]),
      axis.title.y.right = element_text(color = mycolors["trued"]),
      axis.text.y.right = element_text(color = mycolors["trued"]),
      axis.text = element_text(size=24),
      axis.title = element_text(size=26),
      plot.title = element_text(size=34,face = "bold"),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20, face="bold"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1.5),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )

  # if given, plot the smoothed estimates against the reported cases
  if(!is.null(cases.table) && nrow(cases.table) != 0) {
      maxX <- max(cases.table$new_cases_avrg)+0.1*max(cases.table$new_cases_avrg)
      minY <- 0
      maxY <- max(smooth.table$smoothMedian)+0.1*max(smooth.table$smoothMedian)
      ylimMax <- max(cases.table$new_cases_avrg)
      rel_vs_true_ratio <- maxY/maxX
      
      p_smooth_esti_realN <- p_smooth_esti_realN +
      geom_line(aes(x=as.Date(cases.table$date,"%Y-%m-%d"), 
                    y=rescale(cases.table$new_cases_avrg, to=c(minY, maxY), 
                              from=range(cases.table$new_cases_avrg))),size=2.3, color=mycolors["trued"], alpha=0.5) +
        scale_y_continuous(sec.axis = sec_axis(~ . * 1/rel_vs_true_ratio, name = "Reported cases"))
  }
  # #... or only the smoothed line
   else{
     p_smooth_esti_realN <- p_smooth_esti_realN +
       scale_y_continuous(expression(paste(phi[est])))
   }

  if(!is.null(title))
    p_smooth_esti_realN <- p_smooth_esti_realN + labs(title=toString(title))
  
  if(!is.null(outputFile)) {
    ggsave(p_smooth_esti_realN, height = 8, width = 10, dpi = 220, file = outputFile)
  }
  
  #if present, add the point estimates to the plot and save an additional plot
  if(!is.null(point.table) && nrow(point.table) > 0) {
    p_smooth_esti_realN <- p_smooth_esti_realN +
           geom_point(aes(x=as.Date(point.table$meanBinDate,"%Y-%m-%d"), y=point.table$value, size=point.table$sampleSize), colour="dodgerblue4", fill=mycolors["esti"], alpha=0.55) +
           scale_size(name="Sample size",range = c(0,10), breaks = c(100,500,1000,1500,2000)) +
           theme(legend.position="bottom")
  
      if(!is.null(outputFileDots)) {
        ggsave(p_smooth_esti_realN,height = 8, width = 10, dpi = 220, file = outputFileDots)
      }
  }
}
