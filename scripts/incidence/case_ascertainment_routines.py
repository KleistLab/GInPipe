import utils.date_routines as dr
import utils.smoothing_routines as sm

import pandas as pd

def merge_tables(table1, table2):

    if ('date' not in table1.columns) or ('date' not in table2.columns):
        raise ValueError("Error while merging tables. Both tables need a date column.")

    #re-index the date. Get for both tables the same ts
 
    #calculate the date with t=0 (in this case minDate - minDays)
    minDate = min(table1['date'].min(), table2['date'].min())
    date_zero = dr.get_date_zero(1, minDate)
    # filter reported cases information only for dates where a phi value is available
    # note: date in string format is also comparable (lexicographic order = numeric order)
    table1['t'] = [dr.as_days_since_date_zero(d, date_zero) for d in table1['date']]
    table2['t'] = [dr.as_days_since_date_zero(d, date_zero) for d in table2['date']]
    # subset for time where data is available for both data sets
    table1 = table1[table1['t'].isin(table2['t'])]
    table2 = table2[table2['t'].isin(table1['t'])]

    # set for both tables (pandas series) the same index to to rowwise calculations
    table1 = table1.reset_index(drop=True)
    table2 = table2.reset_index(drop=True)

    return pd.merge(table1, table2)


def calculate_minimum_incidence(phi_per_day, rep_cases_per_day):
    
    if len(phi_per_day) != len(rep_cases_per_day):
        raise ValueError("Error while calculating minimum incidence. Smoothed phi and reported cases are not of the same length.")
    
    min_n_true = round(phi_per_day * max(rep_cases_per_day / phi_per_day)).astype(int)
    return min_n_true

