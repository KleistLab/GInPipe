import utils.date_routines as dr
import utils.smoothing_routines as sm

import pandas as pd

def merge_tables(table1, table2):

    if ('date' not in table1.columns) or ('date' not in table2.columns):
        raise ValueError("Error while merging tables. Both tables need a date column.")

    #if ('t' not in table1.columns) or ('t' not in table2.columns):
    #re-index the date. Get for both tables the same ts
    #TODO auslagern 
    #calculate the date with t=0 (in this case minDate - minDays)
    minDate = min(table1['date'].min(), table2['date'].min())
    date_zero = dr.get_date_zero(1, minDate)
    table1['t'] = [dr.as_days_since_date_zero(d, date_zero) for d in table1['date']]
    table2['t'] = [dr.as_days_since_date_zero(d, date_zero) for d in table2['date']]
    # subset for time where data is available for both data sets
    table1 = table1[table1['t'].isin(table2['t'])]
    table2 = table2[table2['t'].isin(table1['t'])]

    # set for both tables (pandas series) the same index to to rowwise calculations
    table1 = table1.reset_index(drop=True)
    table2 = table2.reset_index(drop=True)

    return pd.merge(table1, table2)


def calculate_minimal_incidence(phi_per_day, rep_cases_per_day):
    
    if len(phi_per_day) != len(rep_cases_per_day):
        raise ValueError("Error while calculating minimal incidence. Smoothed phi and reported cases are not of the same length.")
    #calculate maximal scaling coefficient
    #c_t = (rep_cases_per_day / phi_per_day)
    #normalize coefficient
    #c_n = c_t/max(c_t)
    #calculate min_n_true
    #min_n_trueround(rep_cases_per_day/c_n)
    
    min_n_true = round(phi_per_day * max(rep_cases_per_day / phi_per_day)).astype(int)
    return min_n_true


#TODO weg?
def calculate_case_ascertainment_wTests(smoothedPhi_table, testInfo_table, pop = 84334787, sens=1, spec=1):

    #TODO Test: are they pandas.tables? otherwise calculations below won't work
    # only evaluable for the common time
    if 't' not in testInfo_table.columns and not testInfo_table.empty:
        #TODO auslagern? 
        #calculate the date with t=0 (in this case minDate - minDays)
        date_zero = dr.get_date_zero(smoothedPhi_table['t'].iloc[0], smoothedPhi_table['date'].iloc[0])
        testInfo_table['t'] = [dr.as_days_since_date_zero(d, date_zero) for d in testInfo_table['date']] 
    
    testInfo_table = testInfo_table[testInfo_table['t'].isin(smoothedPhi_table['t'])]
    smoothedPhi_table = smoothedPhi_table[smoothedPhi_table['t'].isin(testInfo_table['t'])]

    # set for both tables (pandas series) the same index to to rowwise calculations
    testInfo_table = testInfo_table.reset_index(drop=True)
    smoothedPhi_table = smoothedPhi_table.reset_index(drop=True)
  
    #test positivity
    r_pos = testInfo_table['posRate']
  
    #number of tests
    nt = testInfo_table['newTests']

    phi = smoothedPhi_table['phiSmoothed']

    p_infected_tested = (r_pos - (1-spec))/(sens - (1-spec))
    p_tested = 1-(1-1/pop)**nt
    p_infected = phi/pop

    return pd.DataFrame({'date': smoothedPhi_table['date'], 
        't': smoothedPhi_table['t'],
        'p_infected_tested': p_infected_tested,
        'p_tested': p_tested,
        'p_infected':p_infected,
        'p_tested_infected': p_infected_tested*p_tested/p_infected,
        'pop': pop})


def calculate_minimal_incidence_wTests(caseAscertainment_table, rep_cases_table, 
                                          smoothingWindow=7):

    if 't' not in reportedCases_table.columns:
        #TODO auslagern 
        #calculate the date with t=0 (in this case minDate - minDays)
        date_zero = dr.get_date_zero(caseAscertainment_table['t'][0], caseAscertainment_table['date'][0])
        reportedCases_table['t'] = [dr.as_days_since_date_zero(d, date_zero) for d in reportedCases_table['date']]
    
    # subset for time where data is available for both data sets
    reportedCases_table = reportedCases_table[reportedCases_table['t'].isin(caseAscertainment_table['t'])]
    caseAscertainment_table = caseAscertainment_table[caseAscertainment_table['t'].isin(reportedCases_table['t'])]

    # set for both tables (pandas series) the same index to to rowwise calculations
    reportedCases_table = reportedCases_table.reset_index(drop=True)
    caseAscertainment_table = caseAscertainment_table.reset_index(drop=True)

    #TODO weg? rather rolling average than kernel smoother (the were nans)
    # pti_smoothing = sm.smooth_values(x=caseAscertainment_table['t'], 
    #                             y=caseAscertainment_table['p_tested_infected'], 
    #                             bandwidth=smoothingWindow)
    # case_smoothing = sm.smooth_values(x=caseAscertainment_table['t'], 
    #                             y=caseAscertainment_table['newCases'], 
    #                             bandwidth=smoothingWindow)

    #TODO: direkt mit pandas schöner machen? (like mutate in R)
    caseAscertainment_table['smoothed_P_t_i'] = caseAscertainment_table['p_tested_infected'].rolling(window=smoothingWindow, center=True).mean()
    caseAscertainment_table['newCases'] = reportedCases_table['newCases'].copy(deep=True)
    caseAscertainment_table['smoothedCases'] = caseAscertainment_table['newCases'].rolling(window=smoothingWindow, center=True).mean()

    caseAscertainment_table['pRelDetection'] = caseAscertainment_table['smoothed_P_t_i']/caseAscertainment_table['smoothed_P_t_i'].max()
    caseAscertainment_table['pRelUnderdetection'] = 1.0 - caseAscertainment_table['pRelDetection']

    caseAscertainment_table['minNTrue'] = round(caseAscertainment_table['smoothedCases']/caseAscertainment_table['pRelDetection'])
    caseAscertainment_table['minUnderdetection'] = round(caseAscertainment_table['minNTrue'] - caseAscertainment_table['smoothedCases'])
    # Incidence = per 100.000 per week
    caseAscertainment_table['minTrueIncidence'] = caseAscertainment_table['minNTrue']/caseAscertainment_table['pop']*700000
    caseAscertainment_table['reportedIncidence'] = caseAscertainment_table['smoothedCases']/caseAscertainment_table['pop']*700000

    return caseAscertainment_table