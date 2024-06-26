#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Call smoothing routine for the phi point estimates.
"""

from pathlib import Path

import utils.smoothing_routines as sm
import utils.io_routines as io
import plot.plot_routines as pt
import sys

##################################################
### Read parameters
##################################################

#File containing table with phi estimates per bin
phi_table_file = snakemake.input.phi_table
# Result path
result_file =  snakemake.output[0]
result_path =  str(Path(snakemake.output[0]).parent) 

#suffix for the file name
suffix = ""
if snakemake.params.name:
    suffix = "_" + snakemake.params.name

#Filters:
#minimal bin size
min_bin_size = snakemake.params.min_bin_size
#minimal days to be spanned
min_days_span = snakemake.params.min_days_span
#maximal days to be spanned
max_days_span = snakemake.params.max_days_span
#smoothing window
smoothing_window = snakemake.params.smoothing_bandwidth

try:
    if smoothing_window:
        smoothing_window = int(smoothing_window)
    else:
        smoothing_window = 7
        print("No value given for the smoothing of the phi estimates. Taking default value of " + str(smoothing_window) + ".")

    if min_bin_size:
        min_bin_size =  int(min_bin_size)
    else:
        min_bin_size = 1
        print("No value given for the minimal bin size. Taking default value of " + str(min_bin_size) + ".")


    if min_days_span:
        min_days_span =  int(min_days_span)
    else:
        min_days_span = 1
        print("No value given for the minimal number of days per bin. Taking default value of " + str(min_days_span) + ".")

    if max_days_span:
        max_days_span =  int(max_days_span)
    else:
        # Workaround: if you set max_days_span to empty value (""), 
        # set the max days to at least 3 weeks 
        if max_days_span == "":
            max_days_span = 21
            print("No value given for the maximal number of days per bin. Taking default value of " + str(max_days_span) + ".")
        # If you set it to 0, there is no boundary
        else:
            max_days_span = 0    

except ValueError:
    raise ValueError("Error while smoothing. Filter and smoothing parameters need to be integers.")

##################################################
### Filtering and Smoothing
##################################################

print("-"*80+"\n")
print("Smooth phi estimates with window " + str(smoothing_window)+"\n")
print("-"*80+"\n")

phi_table = io.read_table(phi_table_file)
# filter bins 
phi_table = phi_table[(phi_table["daysPerBin"] >= min_days_span) & \
                               ((not max_days_span) | (phi_table["daysPerBin"] <= max_days_span)) & \
                               (phi_table["sampleSize"] >= min_bin_size)]

if phi_table.empty:
    sys.exit("-"*80+ "\n\n !!!!! No phi estimates to smooth over. Choose different filter criteria. !!!!!\n\n" + "-"*80 + "\n")

smoothed_phi_table = sm.smooth_phi_estimates(phi_table, smoothingBandwidth=smoothing_window)

##################################################
### Plot results
##################################################

out_file = result_path + "/plot_smoothed_phi_estimates" + suffix + ".pdf"
print(" * Plot smoothed phi estimates in " + out_file + "\n")

pt.plot_phi_estimates(phi_table, smoothed_phi_table, out_file)

##################################################
### Write output
##################################################

print(" * Write smoothed phi table to " + result_file + "\n")
io.write_table(smoothed_phi_table, result_file)

