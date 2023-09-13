

from pathlib import Path
import pandas as pd
#from datetime import datetime

import utils.io_routines as io
import utils.smoothing_routines as sm
import incidence.case_ascertainment_routines as ca
import plot.plot_routines as pt
#import sys
# TODO write output and errors to logfile
#with open(snakemake.log[0], "w") as f:
#    sys.stderr = sys.stdout = f


#TODO logging:
# logger = logging.getLogger('ftpuploader')
# hdlr = logging.FileHandler('ftplog.log')
# formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
# hdlr.setFormatter(formatter)
# logger.addHandler(hdlr)
# logger.setLevel(logging.INFO)


print('-'*80 + "\n")
print("*** Calculate the minimal number of infections\n")
print('-'*80 + "\n\n")

##################################################
### Read parameters
##################################################

#File containing table with smoothed phi estimates per day
phi_table_file = snakemake.input.smoothed_phi
#phi_table_file = "/Users/msmith/Documents/RKI/DAKI/GInPipe_HiddenCases/tools/GInPipe/demo_covSonar/results/phi_estimates/smoothed_phi_estimates_testi.csv"
# Result path
mi_result_file =  snakemake.output[0]
result_path =  str(Path(mi_result_file).parent) 

#suffix for the file name
suffix = ""
if snakemake.params.name:
    suffix = "_" + snakemake.params.name

#Reported cases table attributes:
reported_cases_file = snakemake.params.rep_cases
#reported_cases_file = "/Users/msmith/Documents/RKI/DAKI/GInPipe_HiddenCases/tools/ginpipe_covsonar/data/reported_cases_rki_2021.csv"
delimiter = snakemake.params.sep
#","##
col_date = snakemake.params.col_date
#"date"#
col_case = snakemake.params.col_case
#"newCases"#
date_format = snakemake.params.date_format
#"%Y-%m-%d"#

#smoothing window
smoothing_window_mi = snakemake.params.smoothing_bandwidth


#time window (optional)
from_date = snakemake.params.from_date
to_date = snakemake.params.to_date


#TODO
# # name for suffix
# name = snakemake.params.name
# suffix = ""
# if name:
#     suffix = "_"+name


##################################################
### read tables
##################################################

print("*** Read tables")

smoothed_phi_table = io.read_table(phi_table_file)

rep_cases_table = io.read_table(file=reported_cases_file, sep=delimiter)
# Check if the defined column names are given in the reported cases table
if (not col_date in rep_cases_table) or (not col_case in rep_cases_table):
    raise ValueError("Error while reading reported cases. \nColumn " + col_date + 
                " or " + col_case + " does not exist in table " + reported_cases_file + ".\n")

#rename date column of reported cases
rep_cases_table = rep_cases_table.rename(columns={col_date:'date'})



##################################################
### Merge tables and get overlapping time frame
##################################################

print("*** Smooth and filter data points")

# TODO: statt über tabellen gleich ne Klasse machen
mi_table = ca.merge_tables(rep_cases_table, smoothed_phi_table)

##################################################
### Smooth phi estimates and reported cases
##################################################

#TODO smooth phi estimates again? oder eher direkt 21 in the first go??
#TODO so far different naming, stick with smoothed_phi?
mi_table['smoothed_phi'] = sm.smooth_values(x=mi_table['t'],
                                                y=mi_table['phi_smoothed'],
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

print("*** Infer minimal number of infected")
mi_table['min_n_true'] = ca.calculate_minimal_incidence(mi_table['smoothed_phi'], mi_table['smoothed_cases'])

##################################################
### Plot results
##################################################

print("*** Plot results")
out_file = result_path + "/plot_min_infected" + suffix + ".pdf"
pt.plot_min_infected(mi_table, out_file)

##################################################
### Write table
##################################################

print("*** Write result table")
io.write_table(mi_table, file = mi_result_file)

