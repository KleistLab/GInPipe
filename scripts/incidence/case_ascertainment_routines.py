import utils.date_routines as dr
import utils.smoothing_routines as sm

import pandas as pd

class case_ascertainment:
    """
    SAMtoFP class.

    Tranlsate CIGAR strings to SNV strings  
    """

    def __init__(self, rep_cases_table, smoothed_phi_table):
        """
        Calculate minimum number of infected with phi and reported cases.

        :param mi_table: Merged table with phi and reported cases for each day
        :type filename: data.frame
        """
        self.mi_table = merge_tables(rep_cases_table, smoothed_phi_table)

    def smooth_phi(self, bandwidth):
        self.mi_table['smoothed_phi'] = sm.smooth_values(x=self.mi_table['t'],
                                                y=self.mi_table['smoothed_phi'],
                                                bandwidth=bandwidth)['y_smoothed']
        
    def smooth_cases(self, bandwidth):
        self.mi_table['smoothed_cases'] = round(sm.smooth_values(x=self.mi_table['t'],
                                                y=self.mi_table['cases'],
                                                bandwidth=bandwidth)['y_smoothed']).astype(int)
    def calculate_minimum_incidence(self):
        if ('smoothed_cases' not in self.mi_table.columns):
            self.mi_table['smoothed_cases'] = self.mi_table.cases
        self.mi_table['case_phi_ratio'] =  self.mi_table['smoothed_cases']/self.mi_table['smoothed_phi']
        self.mi_table['min_n_true'] = round(self.mi_table['smoothed_phi'] * max(self.mi_table['case_phi_ratio'])).astype(int)

    def get_max_ratio(self):
        if 'case_phi_ratio' not in self.mi_table.columns:
            self.calculate_minimum_incidence()

        return(self.mi_table.iloc[[self.mi_table['case_phi_ratio'].idxmax()]])

        
def merge_tables(table1, table2):

    if ('date' not in table1.columns) or ('date' not in table2.columns):
        raise ValueError("Error while merging tables. Both tables need a date column.")


       # sort tables
    table1 = table1.sort_values(by='date')
    table1 = table1.reset_index(drop=True)
    table2 = table2.sort_values(by='date')
    table2 = table2.reset_index(drop=True)
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

    if len(table1) != len(table2):
        raise ValueError("Error while merging phi and reported cases tables. It seems that a date occurs more than once in one of the tables. Please ensure there is only one entry per date.")

    return pd.merge(table1, table2)

