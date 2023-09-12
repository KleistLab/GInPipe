#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Call smoothing routine for the phi point estimates.
"""

import utils.smoothing_routines as sm
import utils.io_routines as io
#import sys
# TODO write output and errors to logfile
#with open(snakemake.log[0], "w") as f:
#    sys.stderr = sys.stdout = f

##################################################
### Read parameters
##################################################

#File containing table with phi estimates per bin
phi_table_file = snakemake.input.phi_table
# Result path
result_file =  snakemake.output[0]
#Filters:
#minimal bin size
min_bin_size = snakemake.params.min_bin_size
#minimal days to be spanned
min_days_span = snakemake.params.min_days_span
#maximal days to be spanned
max_days_span = snakemake.params.max_days_span
#smoothing window
smoothing_window = snakemake.params.smoothing_window

try:
    if smoothing_window:
        smoothing_window = int(smoothing_window)
    else:
        smoothing_window = 7

    if min_bin_size:
        min_bin_size =  int(min_bin_size)
    else:
        min_bin_size = 1

    if min_days_span:
        min_days_span =  int(min_days_span)
    else:
        min_days_span = 1
    
    if max_days_span:
        max_days_span =  int(max_days_span)
    else:
        max_days_span = 0
except ValueError:
    raise ValueError("Error while smoothing. Filter and smoothing parameters need to be integers.")

##################################################
### Filtering and Smoothing
##################################################

print("*"*80+"\n")
print("*** Smooth phi estimates with window " + str(smoothing_window)+"\n")
print("*"*80+"\n")
phi_table = io.read_table(phi_table_file)
# filter bins 
phi_table = phi_table[(phi_table["daysPerBin"] >= min_days_span) & \
                               ((not max_days_span) | (phi_table["daysPerBin"] <= max_days_span)) & \
                               (phi_table["sampleSize"] >= min_bin_size)]


smoothed_phi_table = sm.smooth_phi_estimates(phi_table, smoothingBandwidth=smoothing_window)

##################################################
### Write output
##################################################

print("\n...write smoothed phi table to " + result_file + "\n")
io.write_table(smoothed_phi_table, result_file)

#TODO plot
