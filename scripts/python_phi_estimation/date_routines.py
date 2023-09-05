
from datetime import datetime, timedelta
from functools import singledispatch


#%j 	Day of the year
#%u Weekday as a decimal number (1–7, Monday is 1).
#%U Week of the year as decimal number (00–53) using Sunday as the first day 1 of the week (and typically with the first Sunday of the year as day 1 of week 1). The US convention.
#%W Week 00-53 with Monday as first day of the week


#constants
DATE_FORMAT = '%Y-%m-%d'

#' Add the number of days to the date
#'
#' @param days as integer
#' @param date_str string representation of date in format %Y-%m-%d
#' @return date as string
#' @export
@singledispatch
def add_days_to_date(days, date_str):
	date = datetime.strptime(date_str, DATE_FORMAT) + timedelta(days=int(days))
	return date.strftime(DATE_FORMAT)

# @add_days_to_date.register
# def _(days:int, date_str:str):
# 	date = datetime.strptime(date_str, DATE_FORMAT) + timedelta(days)
# 	return date.strftime(DATE_FORMAT)

# #TODO WORKAROUND because union typing is not working
# add_days_to_date.register
# def _(days:float, date_str:str):
# 	date = datetime.strptime(date_str, DATE_FORMAT) + timedelta(days)
# 	return date.strftime(DATE_FORMAT)

@add_days_to_date.register
def _(days:list, date_str:str):
	return [add_days_to_date(d, date_str) for d in days]

#' Substract day index t from date to get the date for day zero.
#'
#' @param t days since day zero 
#' @param date_str date of the "t"th day
#' @return date at day zero
#' @export
def get_date_zero(t: int|float, date_str: str):
	return add_days_to_date(-t, date_str)
	

# TODO: noch checken ob einzeldatum oder vektor... 
#' Return the days since the date of day zero for the given date.
#'
#' @param date 
#' @param date zero
#' @return time delta in days (integer)
#' @export
@singledispatch
def as_days_since_date_zero(date, date_zero):
    ... #else case TODO: error handling

@as_days_since_date_zero.register
def _(date:str, date_zero:str):
	time_delta = datetime.strptime(date, DATE_FORMAT) - datetime.strptime(date_zero, DATE_FORMAT)
	return time_delta.days

@as_days_since_date_zero.register
def _(date:list, date_zero:str):
   return [as_days_since_date_zero(d, date_zero) for d in date]



#' Return the date in string format as date in datetime format.
@singledispatch
def str_to_date(date_str):
    ... #else case TODO: error handling

@str_to_date.register
def _(date_str:str):
    return datetime.strptime(date_str, DATE_FORMAT)

@str_to_date.register
def _(date_str:list):
   return [str_to_date(d) for d in date_str] 



# #' Return the calender week of the date (ISO 8601)
# #'
# #' @param date to be transformed
# #' @return the numeric calendar week
# #' @export
# def as_cal_week(date):
#   return(as.numeric(format(date, "%W")))



# #' Return the date from the day of the year of the date as yyyyddd.
# #'
# #' @param doy day of the year in the format "yyyyddd"
# #' @return the date
# #' @export
# doyAsDate <- function(doy) {
#   return(as.Date(x = as.character(doy), format = "%Y%j"))
# }

# #' Return for each date in the vector the days since the given date "minDate".
# #'
# #' @param dates a vector of dates
# #' @param minDate a date
# #' @return a vector of days (integers)
# #' @export
# asDays_since_global_d0 <- function(dates, minDate) {
#   dates <- as.Date(dates)
#   #minDate <- min(dates)
#   return(as.numeric(dates - minDate))
# }

# #' Return the dates when adding the given vector of integer days to the given date "minDate"
# #'
# #' @param days a vector of days (integers)
# #' @param minDate a date
# #' @return a vector of dates
# #' @export
# daysAsDate <- function(days, minDate) {
#   return(as.Date(as.Date(minDate) + days))
# }

# #' Return a representative date of the calender week. Here: the last day.
# #'
# #' @param week a calendar week (integer)
# #' @param year a year (integer) in the format yyyy
# #' @return the date of the last day of the calendar week
# #' @export
# asLastDayOfWeek <- function(week, year) {
#   return(as.Date(x = paste0(year, "-",week,"-7"), format = "%Y-%W-%u"))
# }
