#%j 	Day of the year
#%u Weekday as a decimal number (1–7, Monday is 1).
#%U Week of the year as decimal number (00–53) using Sunday as the first day 1 of the week (and typically with the first Sunday of the year as day 1 of week 1). The US convention.
#%W Week 00-53 with Monday as first day of the week


#' Return the calender week of the date (ISO 8601)
#'
#' @param date to be transformed 
#' @return the numeric calendar week
as.calWeek <- function(date) {
  return(as.numeric(format(date, "%W")))
}


#' Return the day of the year of the date as "yyyyddd".
#'
#' @param date to be transformed
#' @return the numeric day of the year in the format "yyyyddd"
as.dayOfYear <- function(date) {

  return(as.numeric(format(date, "%Y%j")))
}

#' Return the date from the day of the year of the date as yyyyddd.
#'
#' @param doy day of the year in the format "yyyyddd"
#' @return the date 
doy.as.Date <- function(doy) {
  return(as.Date(x = as.character(doy), format = "%Y%j"))
}

#' Return for each date the days since the minimum date in the vector.
#'
#' @param dates a vector of dates
#' @return a vector of days (integers)
as.days_since_d0 <- function(dates) {
  dates <- as.Date(dates)
  minDate <- min(dates)
  return(as.numeric(dates - minDate))
}

#' Return for each date in the vector the days since the given date "minDate".
#'
#' @param dates a vector of dates
#' @param minDate a date
#' @return a vector of days (integers)
as.days_since_global_d0 <- function(dates, minDate) {
  dates <- as.Date(dates)
  #minDate <- min(dates)
  return(as.numeric(dates - minDate))
}

#' Return the dates when adding the given vector of integer days to the given date "minDate"
#'
#' @param days a vector of days (integers)
#' @param minDate a date
#' @return a vector of dates 
days.as.Date <- function(days, minDate) {
  return(as.Date(as.Date(minDate) + days))
}

#' Return a representative date of the calender week. Here: the last day.
#'
#' @param week a calendar week (integer)
#' @param year a year (integer) in the format yyyy
#' @return the date of the last day of the calendar week
as.lastDayOfWeek <- function(week, year) {
  return(as.Date(x = paste0(year, "-",week,"-7"), format = "%Y-%W-%u"))
}
