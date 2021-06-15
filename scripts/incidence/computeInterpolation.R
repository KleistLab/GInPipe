#!/usr/local/bin Rscript --vanilla --slave

# clear the environment variables
rm(list = ls())

options(warn=-1)

# install needed packages
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")"))))
    return(TRUE)
  install.packages(package,repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

#"ggformula"
for(p in c("ggplot2","MASS","R0","scales")) {
  dynamic_require(p)
}

#Read in arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  cat("\nCall the script with 3 base arguments and 5 optional arguments: estimatesFile label outputFile R0 (reportedCasesFile delim dateColumName newCasesColumnName dateFormat )\n\n
1. The estimates file contains a tab separated table with headers and has at least 3 columns: \n
t value variance.\n
2. A label for the data set, e.g. a country or city.\n
3. The output path. The output tables are written to the given directory,
      which is created if it does not exist yet.\n\n
Optionally: \n
4. The reported cases file contains a table with reported cases on each date. \n
The separator and the column names can be chosen arbitrarily and are defined with the following parameters.\n
5. Delim gives the delimiter in the reported cases table. \n
6. The column name for the the date in the reported cases table. \n
7. The column name for the the number of cases in the reported cases table. \n
8. The date format in table.\n\n")
  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
inputFile<-normalizePath(args[1])
country = toString(args[2])
outputFile<-file.path(args[3])
table_name = NA
table_delim = NA
table_date_col = NA
table_active_col = NA
table_date_format = NA
if(length(args)==8) {
  table_name = args[4]
  table_delim = args[5]
  table_date_col = args[6]
  table_active_col = args[7]
  table_date_format = args[8]
}


# A pre-set date format for metric output
input_date_format = "%Y-%m-%d"

# Read input
input.table = read.table(inputFile, header=T, sep = "\t")
# Define output file and output directory
fileName<-basename(outputFile)
outputDir<- dirname(outputFile)
if(!dir.exists(outputDir))
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
outputDir <- normalizePath(outputDir)
outputFile<-paste0(outputDir,"/",fileName)


# Read reporte cases, if existing
cases.table.full = data.frame(matrix(ncol=0,nrow=0))
cases.table = data.frame(matrix(ncol=0,nrow=0))

# Some fake date
minDate1 <- integer(0)
class(minDate1) <- "Date"

if (!is.na(table_name) && file.exists(table_name)) {
  cases.table.full = read.table(table_name, header=T, sep=table_delim)
  # change column names as they can be chosen flexibly, but we require date and new_cases
  names(cases.table.full)[names(cases.table.full) == table_date_col] <- "date"
  names(cases.table.full)[names(cases.table.full) == table_active_col] <- "new_cases"
  minDate1 = min(as.Date(cases.table.full$date, table_date_format))
  #minDates <- c(minDates,minDate1)
  cases.table.full$date <- as.Date(cases.table.full$date, table_date_format)
  # Make cases table with t, new_cases, dates
  cases.table <- data.frame(new_cases=cases.table.full$new_cases,date=cases.table.full$date)
  if (!"new_cases" %in% colnames(cases.table)) {
    cat("\n Reported cases table has no column named 'new_cases'. Rename existing columns
                or provide a new table\n")
    #terminate without saving workspace
    quit("no")}
  cases.table$new_cases[cases.table$new_cases<0] <- 0
  cases.table$new_cases[is.na(cases.table$new_cases)] <- 0
  cases.table$new_cases_avrg <- filter(cases.table$new_cases, rep(1/7,7))
  cases.table <- na.omit(cases.table)
  cases.table$country <- rep(country,nrow(cases.table))
}


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
this.dir <- getScriptPath()

#getSrcDirectory(function(x) {x})
#if (rstudioapi::isAvailable()) {
#  this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#}else {
#  this.dir <- getScriptPath()
#}
setwd(this.dir)
# load other r routines

source("dateFormatRoutines.R")
source("splineRoutines.R")
source("plotRoutines.R")
# minDate to use in all plotted tables
minDate2 = min(as.Date(input.table$meanBinDate,input_date_format))
minDates <- c(minDate1,minDate2)
minDate = min(minDates)
# Compute the splines and dot sizes
cat("--- Compute interpolation ---\n\n")
input.table$t <- as.days_since_global_d0(input.table$meanBinDate,minDate)
if (!is.na(table_name)) {cases.table$t <- as.days_since_global_d0(cases.table$date,minDate)}
pointSize <- c()
for (i in (1:nrow(input.table))) {
  # rand*(UB1-LB1) + LB1
  point <- ((10-1)*(input.table$sampleSize[i]-
            min(input.table$sampleSize)))/(max(input.table$sampleSize)-
            min(input.table$sampleSize)) + 1
  pointSize <- c(pointSize, point)
}
input.table$pointSize <- pointSize

# Replace negative number of new cases with 0 - happens if calculated from cumulative confirmed
# cases count
interp.table <-  computeInterpolation(input.table, seq(min(input.table$t), max(input.table$t)), input.table$sampleSize)
# Plot interpolated curve with 95% CI starting on global minDate
interp.table["date"] = days.as.Date(interp.table$t, minDate)
#interp.table <- interp.table[!is.na(interp.table$median) ,]
interp.table[is.na(interp.table)] = 0
# Remove rows with zeros in smooth median
interp.table <- interp.table[interp.table$smoothMedian != 0,]
# Add date column
interp.table$date <- days.as.Date(interp.table$t, minDate)
# Write tables and plot
write.csv(interp.table,paste0(outputDir,"/",fileName), row.names = F, col.names = T)
outputFileInter<-paste0(outputDir,"/",substr(fileName,1,(nchar(fileName)-4)),".pdf")
outputFileInterDots<-paste0(outputDir,"/","wdots_",substr(fileName,1,(nchar(fileName)-4)),".pdf")
plotInterpolationWithNewCases(cases.table, interp.table, input.table, minDate, outputFileInter, outputFileInterDots,country)

# Write estimated theta table
input.table.clean <- data.frame(date=input.table$meanBinDate,phi=input.table$value,variance=input.table$variance)
write.csv(input.table.clean,paste0(outputDir,"/phi_estimate_",country,".csv"), row.names = F)
