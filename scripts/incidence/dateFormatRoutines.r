#%j 	Day of the year
#%u Weekday as a decimal number (1–7, Monday is 1).
#%U Week of the year as decimal number (00–53) using Sunday as the first day 1 of the week (and typically with the first Sunday of the year as day 1 of week 1). The US convention.
#%W Week 00-53 with Monday as first day of the week

# return the calender week of the date (ISO 8601)
as.calWeek <- function(date) {
  return(as.numeric(format(date, "%W")))
}


# return the day of the year of the date
as.dayOfYear <- function(date) {

  return(as.numeric(format(date, "%Y%j")))
}

doy.as.Date <- function(doy) {
  return(as.Date(x = as.character(doy), format = "%Y%j"))
}

as.days_since_d0 <- function(dates) {
  dates <- as.Date(dates)
  minDate <- min(dates)
  return(as.numeric(dates - minDate))
}

as.days_since_global_d0 <- function(dates,minDate) {
  dates <- as.Date(dates)
  #minDate <- min(dates)
  return(as.numeric(dates - minDate))
}

days.as.Date <- function(days, minDate) {
  return(as.Date(as.Date(minDate) + days))
}

# return a representative of the calender week. For now: the last day
as.lastDayOfWeek <- function(week, year) {
  return(as.Date(x = paste0(year, "-",week,"-7"), format = "%Y-%W-%u"))
}
#
# fix.date.format <- function(dates) {
#   dates.str <- c()
#   date.format <- "%Y-%m-%d"
#   date.trunc1 <- "%Y-%m"
#   date.trunc2 <- "%Y"
#   for (i in (1:length(dates))) {
#     if (is.na(as.Date(dates[i], date.format))) {
#       print(dates[i])
#       if (!is.na(as.Date(dates[i], date.trunc1))) {
#         date <- as.Date(paste0(toString(dates[i]), "-15"), format = "%Y-%m-%d")
#         print(date)
#         dates.str <- c(dates.str,date)
#       }
#       else if (!is.na(as.Date(dates[i], date.trunc2))) {
#         date <- as.Date(paste0(toString(dates[i]), "-03-01"), format = "%Y-%m-%d")
#         print(date)
#         dates.str <- c(dates.str,date)
#       }
#     }
#     else {
#       dates.str <- c(dates.str, toString(dates[i]))
#     }
#   }
#   return(dates)
# }
