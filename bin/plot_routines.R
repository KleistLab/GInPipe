#!/usr/local/bin Rscript --vanilla --slave

# install needed packages
dynamic_require <- function(package){
    if(eval(parse(text=paste("require(",package,")")))) {
        return(TRUE)
    }
    install.packages(package)
    return(eval(parse(text=paste("require(",package,")"))))
}

for(p in c("ggplot2", "scales", "devtools")) {
    dynamic_require(p)
}

library(ggplot2)
library(scales)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

phi_file <- args[1]
trueN_file <- args[2]

phi.tab <- read.table(phi_file, header=T, sep='\t')
trueN.tab <- read.table(trueN_file, header=T, sep='\t')

#phi.est <- as.integer(phi.tab$phi)
#trueN <- as.numeric(trueN.tab$true_N)
print(phi.tab)
print(trueN.tab)

plot_sim_result <- function(phi,trueN,out) {
    mycolors <- c("estid"="dodgerblue4", "trued"="darkred","esti"="dodgerblue4", "true"="red") # nolint
    p <- ggplot() +
        geom_line(aes(x = seq(1,nrow(phi)), y = phi$phi), colour=mycolors["esti"]) + # nolint
        geom_line(aes(x = seq(1,nrow(trueN)), y = trueN$true_N), colour=mycolors["true"]) +  # nolint
        xlab("t") +  # nolint
        ylab("N") +
        theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # nolint
          panel.background = element_blank(), axis.line = element_line(colour = "black")) # nolint

    ggsave(p,height = 8, width = 10, dpi = 220, file = out)
}

outfile = "outplot.pdf"
plot_sim_result(phi.tab,trueN.tab,outfile)