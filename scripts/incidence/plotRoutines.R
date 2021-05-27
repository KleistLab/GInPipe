#library('latex2exp')
library(ggplot2)
library(grid)
source("dateFormatRoutines.r")

plotR0Package <- function(rzero.table,method,outputFile){
  p_interp_d1 <- ggplot() +
    geom_line(aes(x=as.Date(rzero.table$date,"%Y-%m-%d"), y=rzero.table$value), colour="darkred", size=2, alpha=0.7)+
    geom_line(aes(x=as.Date(rzero.table$date,"%Y-%m-%d"), y=rzero.table$lower), size=0.5, alpha=0.7, linetype="dashed")+
    geom_line(aes(x=as.Date(rzero.table$date,"%Y-%m-%d"), y=rzero.table$upper), size=0.5, alpha=0.7, linetype="dashed")+
    #geom_point(aes(tt.df[,1], d1_gam), alpha=0.5, size=2) +
    #ylim(c(0,100))+
    xlab("") +
    scale_x_date(date_breaks = "2 months" , date_labels = "%b %Y") +
    ylab(expression(paste(R[0]))) +
    ggtitle(method) +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
    ggsave(p_interp_d1,
             height = 6,
             width = 10,
             dpi = 220,

             file = outputFile)
}

plotSplineWithNewCases <-function(data.table, input.table, spline.table, outputFile, minDate) {
  minX <- max(min(data.table$new_cases),0)
  maxX <- max(data.table$new_cases)
  minY <- max(min(spline.table$value),0)
  maxY <- max(spline.table$value)
  ylimMax <- max(data.table$new_cases)
  mycolors <- c("estid"="dodgerblue4", "trued"="darkred","esti"="dodgerblue4", "true"="red")
  p_spline_esti_realN <- ggplot() +
    #geom_point(aes(doy.as.Date(data.table$doy), data.table$trueN*rel_vs_true_ratio), size=2, color=mycolors["trued"], alpha=0.5)+
    geom_line(aes(days.as.Date(data.table$t, minDate), rescale(data.table$new_cases, minX, maxX, minY, maxY)),size=2, color=mycolors["trued"], alpha=0.5)+
    geom_point(aes(x=days.as.Date(input.table$t, minDate), y=input.table$value), colour=mycolors["esti"], size=input.table$pointSize, alpha=0.3)+
    #geom_line(aes(x=doy.as.Date(spline.table$doy), y=spline.table$value_trueN*rel_vs_true_ratio), color=mycolors["true"], size = 2, alpha=0.5)+
    #geom_line(aes(x=days.as.Date(spline.table$t, minDate), y=rescale(spline.table$value_trueN, minX, maxX, minY, maxY)), color=mycolors["true"], size = 2, alpha=0.5)+
    geom_line(aes(x=days.as.Date(spline.table$t, minDate), y=spline.table$value), color=mycolors["esti"], size = 2, alpha=0.5)+
    scale_size_continuous(range = c(min(input.table$pointSize), max(input.table$pointSize))) +
    scale_y_continuous(
      expression(paste(theta[est])),
      #sec.axis = sec_axis(~ . * 1/rel_vs_true_ratio, name = "new cases")
      sec.axis = sec_axis(~ rescale(., minY, maxY, minX, maxX ), name = "New cases")) +
    xlab("")+
    #coord_cartesian(xlim=c(minDate, maxDate),ylim=c(0,ylimMax)) +
    scale_x_date(date_breaks = "2 months" , date_labels = "%b %Y") +
    theme(
      axis.text.x=element_text(angle = 0, vjust = 1, hjust=0),
      axis.title.y = element_text(color = mycolors["estid"]),
      axis.text.y = element_text(color = mycolors["estid"]),
      axis.title.y.right = element_text(color = mycolors["trued"]),
      axis.text.y.right = element_text(color = mycolors["trued"]),
      axis.text = element_text(size=20),
      axis.title = element_text(size=20),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20, face="bold"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )

  ggsave(p_spline_esti_realN,
           height = 6,
           width = 10,
           dpi = 220,

           file = outputFile)
}


plotInterpolationWithNewCases <- function(cases.table,interp.table,input.tables,minDate,outputFile,outputFileDots,group) {
  minX <- 0
  maxX <- max(cases.table$new_cases_avrg)+0.1*max(cases.table$new_cases_avrg)
  minY <- 0
  maxY <- max(interp.table$smoothMedian)+0.1*max(interp.table$smoothMedian)
  ylimMax <- max(cases.table$new_cases_avrg)
  rel_vs_true_ratio = maxY/maxX
  country = toString(group)
  print(country)
  # Min/Max dates
  #xlimMax = min(max(days.as.Date(interp.table$t, minDate)),max(days.as.Date(cases.table$t, minDate)))
  #xlimMin = max(min(days.as.Date(interp.table$t, minDate)),min(days.as.Date(cases.table$t, minDate)))
  mycolors <- c("estid"="dodgerblue4", "trued"="darkred","esti"="dodgerblue4", "true"="red")
  p_spline_esti_realN <- ggplot()
    #rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE), ...)
    #geom_histogram(data=meta.table, aes(x=days.as.Date(meta.table$t, minDate),y=..density..), fill="black", alpha=0.2, bins=round(nrow(input.table)/7)) +
    #geom_point(aes(doy.as.Date(data.table$doy), data.table$trueN*rel_vs_true_ratio), size=2, color=mycolors["trued"], alpha=0.5)+
    if(nrow(cases.table) != 0) {
      p_spline_esti_realN <- p_spline_esti_realN +
      geom_line(aes(x=days.as.Date(cases.table$t, minDate), y=rescale(cases.table$new_cases_avrg, to=c(minY, maxY), from=range(cases.table$new_cases_avrg))),size=2.3, color=mycolors["trued"], alpha=0.5)+
      geom_line(aes(x=days.as.Date(interp.table$t, minDate), y=interp.table$smoothMedian), size=2.3, colour=mycolors["esti"], alpha=0.5)+
      geom_line(aes(x=days.as.Date(interp.table$t, minDate), y=interp.table$smooth5), size=1.3, colour=mycolors["esti"], alpha=0.55, linetype="dashed")+
      #geom_rug(data=meta.table, aes(x = days.as.Date(meta.table$t, minDate)), inherit.aes = F, alpha=0.05)+
      #geom_point(aes(x = days.as.Date(meta.table$t, minDate), y = 0),colour="black",size=3,shape=1,alpha=0.1) +
      geom_line(aes(x=days.as.Date(interp.table$t, minDate), y=interp.table$smooth95), size=1.3, colour=mycolors["esti"], alpha=0.55, linetype="dashed")+
      scale_y_continuous(
        expression(paste(phi[est])),
        sec.axis = sec_axis(~ . * 1/rel_vs_true_ratio, name = "Reported cases")) +
        #sec.axis = sec_axis(~ rescale(cases.table$new_cases_avrg, to=c(minY, maxY), from=range(cases.table$new_cases_avrg)), name = "Reported cases")) +
      xlab("")+
      labs(title=country)+
      #coord_cartesian(xlim=c(minDate, maxDate),ylim=c(0,ylimMax)) +
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
    }
    else{
      p_spline_esti_realN <- p_spline_esti_realN +
      geom_line(aes(x=days.as.Date(interp.table$t, minDate), y=interp.table$smoothMedian), size=2.3, colour=mycolors["esti"], alpha=0.5)+
      geom_line(aes(x=days.as.Date(interp.table$t, minDate), y=interp.table$smooth5), size=1.3, colour=mycolors["esti"], alpha=0.55, linetype="dashed")+
      #geom_rug(data=meta.table, aes(x = days.as.Date(meta.table$t, minDate)), inherit.aes = F, alpha=0.05)+
      #geom_point(aes(x = days.as.Date(meta.table$t, minDate), y = 0),colour="black",size=3,shape=1,alpha=0.1) +
      geom_line(aes(x=days.as.Date(interp.table$t, minDate), y=interp.table$smooth95), size=1.3, colour=mycolors["esti"], alpha=0.55, linetype="dashed")+
      scale_y_continuous(expression(paste(phi[est]))) +
      xlab("")+
      labs(title=country)+
      #coord_cartesian(xlim=c(minDate, maxDate),ylim=c(0,ylimMax)) +
      scale_x_date(date_breaks = "1 months", date_labels = "%b %Y") +
      theme(
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(color = mycolors["estid"]),
        axis.text.y = element_text(color = mycolors["estid"]),
        axis.text = element_text(size=24),
        axis.title = element_text(size=26),
        plot.title = element_text(size=34,face = "bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
      )
    }
    ggsave(p_spline_esti_realN,
             height = 8,
             width = 10,
             dpi = 220,
             file = outputFile)
    p_spline_esti_realN <- p_spline_esti_realN +
         geom_point(aes(x=days.as.Date(input.table$t, minDate), y=input.table$value, size=input.table$sampleSize), colour="dodgerblue4", fill=mycolors["esti"], alpha=0.55) +
         scale_size(name="Sample size",range = c(0,10), breaks = c(100,500,1000,1500,2000)) +
         theme(
           axis.text.x=element_text(size = 20, angle = 45, vjust = 1, hjust=1),
           axis.title.y = element_text(color = mycolors["estid"]),
           axis.text.y = element_text(color = mycolors["estid"]),
           axis.title.y.right = element_text(color = mycolors["trued"]),
           axis.text.y.right = element_text(color = mycolors["trued"]),
           axis.text = element_text(size=24),
           axis.title = element_text(size=26),
           plot.title = element_text(size=34,face = "bold"),
           legend.position="bottom",
           legend.text = element_text(size=20),
           legend.title = element_text(size=20, face="bold"),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "black", fill=NA, size=1.5),
           panel.background = element_blank(), axis.line = element_line(colour = "black")
         )

    ggsave(p_spline_esti_realN,
             height = 8,
             width = 10,
             dpi = 220,
             file = outputFileDots)
  }

plotSpline <- function(input.table, spline.table, outputFile) {
  minDate <- min(as.Date(input.table$meanBinDate))
  p_spline <- ggplot() +
    geom_ribbon(aes(x=days.as.Date(spline.table$t, minDate), ymin = spline.table$lowerSim, ymax = spline.table$upperSim), fill = "grey70",alpha=0.5) +
    geom_point(aes(days.as.Date(input.table$t, minDate), y=input.table$value, size=1/input.table$variance), alpha=0.5, show.legend = F) +
    geom_line(aes(x=days.as.Date(spline.table$t,minDate), spline.table$value), color="red", size = 2, alpha=0.7) +
    xlab("") +
    scale_x_date(date_breaks = "months" , date_labels = "%Y-%m-%d") +
    ylab(expression(paste("estim. ", theta)))+
    theme(axis.text = element_text(size=9,face="bold"),
          axis.title = element_text(size=12),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))


  ggsave(p_spline,
       height = 8,
       width = 16,
       dpi = 220,
       device = "pdf",
       file = outputFile)
}
