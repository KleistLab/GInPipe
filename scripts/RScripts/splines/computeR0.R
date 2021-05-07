#!/usr/local/bin Rscript --vanilla --slave

# clear the environment variables
rm(list = ls())


# install needed packages
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")"))))
    return(TRUE)
  install.packages(package,repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

#"ggformula"
for(p in c("ggplot2","R0")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  cat("\nCall the script with two arguments: \n\n
    1. Interpolation table. The table should contain columns t, date, smoothMedian; \n
    2. Output file path")
  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-args[1]
outputFile<-args[2]
# read input
input.table = read.table(inputFile, header=T, sep = ",")
# Output file
fileName<-basename(outputFile)
outputDir<- dirname(outputFile)
print(outputFile)
if(!dir.exists(outputDir))
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
outputDir <- normalizePath(outputDir)
outputFile<-paste0(normalizePath(outputDir),"/",fileName)

### All absolute paths are set, change working directory
# set working directory to call other R Scripts
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
this.dir <- getSrcDirectory(function(x) {x})
if (rstudioapi::isAvailable()) {
  this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}else {
  this.dir <- getScriptPath()
}
setwd(this.dir)
# load other r routines
source("plotRoutines.r")

cat("--- Compute R0 ---\n\n")
# R0 package
# Wallinga and Teunis (2004)
# Generation intervals distribution - gamma with mean = 5, sd = 1
GT <- generation.time(type = "gamma",
                      val = c(5,1), truncate = NULL, step = 1, first.half = TRUE,
                      p0 = TRUE)
# Pseudocount
input.table$smoothMedian <- input.table$smoothMedian+1
# Compute R0
td <- est.R0.TD(as.numeric(unlist(round(input.table$smoothMedian))),GT=GT,t=input.table$date)
td.table <- data.frame(t=input.table[td$begin.nb:td$end.nb,]$t,value=as.vector(td$R),lower=as.vector(td$conf.int$lower),upper=as.vector(td$conf.int$upper))
# plot
td.table$date <- input.table[td$begin.nb:td$end.nb,]$date
outputFileWT <- paste0(normalizePath(outputDir),"/",substr(fileName,1,(nchar(fileName)-4)),".pdf")
plotR0Package(td.table,"WT04",outputFileWT)
# Write table
write.csv(td.table,paste0(outputDir,"/r0.csv"), row.names = F, col.names = T)
