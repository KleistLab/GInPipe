#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert BAM file into dna profile table and write table.
"""

from pathlib import Path
import pandas as pd

import utils.io_routines as io
import utils.smoothing_routines as sm
import incidence.case_ascertainment_routines as ca
import plot.plot_routines as pt


print('-'*80 + "\n")
print("Calculate the minimal number of infections\n")
print('-'*80 + "\n\n")

##################################################
### Read parameters
##################################################

#File containing table with smoothed phi estimates per day
phi_table_file = snakemake.input.smoothed_phi
# Result path
mi_result_file =  snakemake.output[0]
result_path =  str(Path(mi_result_file).parent) 

#suffix for the file name
suffix = ""
if snakemake.params.name:
    suffix = "_" + snakemake.params.name

#Reported cases table attributes:
reported_cases_file = snakemake.params.rep_cases
delimiter = snakemake.params.sep
col_date = snakemake.params.col_date
col_case = snakemake.params.col_case
date_format = snakemake.params.date_format

#smoothing window
smoothing_window_mi = snakemake.params.smoothing_bandwidth


#time window (optional)
from_date = snakemake.params.from_date
to_date = snakemake.params.to_date

##################################################
### read tables
##################################################

print(" * Read tables")

smoothed_phi_table = io.read_table(phi_table_file)

rep_cases_table = io.read_table(file=reported_cases_file, sep=delimiter)
# Check if the defined column names are given in the reported cases table
if (not col_date in rep_cases_table) or (not col_case in rep_cases_table):
    raise ValueError("Error while reading reported cases. \nColumn " + col_date + 
                " or " + col_case + " does not exist in table " + reported_cases_file + ".\n")

#rename date column of reported cases
rep_cases_table = rep_cases_table.rename(columns={col_date:'date'})

# sort tables
smoothed_phi_table = smoothed_phi_table.sort_values(by='date')
smoothed_phi_table = smoothed_phi_table.reset_index(drop=True)
rep_cases_table = rep_cases_table.sort_values(by='date')
rep_cases_table = rep_cases_table.reset_index(drop=True)

##################################################
### Merge tables and get overlapping time frame
##################################################

print(" * Smooth and filter data points")

# TODO: Refactoring: statt über tabellen gleich ne Klasse machen
mi_table = ca.merge_tables(rep_cases_table, smoothed_phi_table)
##################################################
### Smooth phi estimates and reported cases
##################################################

mi_table['smoothed_phi'] = sm.smooth_values(x=mi_table['t'],
                                                y=mi_table['smoothed_phi'],
                                                bandwidth=smoothing_window_mi)['y_smoothed']

mi_table['smoothed_cases'] = sm.smooth_values(x=mi_table['t'],
                                                y=mi_table[col_case],
                                                bandwidth=smoothing_window_mi)['y_smoothed']

##################################################
### Filter by dates
##################################################

if(from_date):
    mi_table = mi_table.loc[(mi_table['date'] >= from_date)]
if(to_date):
    mi_table = mi_table.loc[(mi_table['date'] <= to_date)]
               
##################################################
### Calculate minimal true infected
##################################################

print(" * Infer minimal number of infected")
mi_table['min_n_true'] = ca.calculate_minimal_incidence(mi_table['smoothed_phi'], mi_table['smoothed_cases'])

##################################################
### Plot results
##################################################

print(" * Plot results")
out_file = result_path + "/plot_min_infected" + suffix + ".pdf"
pt.plot_min_infected(mi_table, out_file)

##################################################
### Write table
##################################################

print(" * Write result table")
io.write_table(mi_table, file = mi_result_file)

