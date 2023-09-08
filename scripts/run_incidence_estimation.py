#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main call to run binning and phi estimation. 
Writing a table containing the estimated phi value for each bin and sequence statistics per days.
TODO fertig beschreiben/anpassen
#' @param snv_file File containing a table in covSonar-output.... 
#' I.e. the sequences are represented as string of single nucleotide variants (SNVs) seperated by a blank, 
#' and a date is given in the format yyyy-mm-dd.
#' @param freq_threshold Threhsold for filtering. Positions which don't occur more often than the threshold in the given sequence set are filtered out.
#TODO @param filter_mutation (optional) Boolean indicating if the threshold should be applied to each particular SNV instead of the position. (default = FALSE)
#' @param masked_positions (optional) Vector containing positions which should be masked/filtered out for the genotyping.
#' @param days_per_bin Creating consecutive sequence bins spanning the given number of days.
#' @param seqs_per_bin Creating consecutive sequence bins containing the given number of sequences.
#' @param result_path(optional)  
#' 
"""

import os
from pathlib import Path
import numpy as np
#import csv
import pandas as pd
#import datetime
#import random

import utils.io_routines as io
import phi.snv_filter_routines as filt
import phi.binning_routines as bin
import phi.seq_info_routines as si
import sys

# write output and errors to logfile
#with open(snakemake.log[0], "w") as f:
#    sys.stderr = sys.stdout = f


# Output path
result_path =  str(Path(snakemake.output[0]).parent) 
#result_path=Path("/Users/msmith/Documents/RKI/DAKI/test_results_ginSonar")

#suffix for the file name
suffix = ""
if snakemake.params.name:
    suffix = "_" + snakemake.params.name

#File containing table with Single nucleotide variants (snv) per sequence and date
snv_file = snakemake.input.snv_file
#snv_file = "/Users/msmith/Documents/RKI/DAKI/GInPipe_HiddenCases/tools/ginpipe_covsonar/data/covSonar_accessions_Germany_2021.csv"

# Grouping variable - suffix for all files
#file_suffix = snakemake.params.group

# Base frequency cutoff
freq_cutoff = snakemake.params.cutoff
#freq_cutoff = 2

# Masking parameters
masking_file = snakemake.params.vcf_file
masked_positions_str = snakemake.params.masked_positions
#masked_positions_str = "1-10,29700-30000"


# Binning parameters
days_per_bin = snakemake.params.days_per_bin
seq_per_bin = snakemake.params.seq_per_bin
#days_per_bin = [7,14]
#seq_per_bin = [1000, 2000, 5000]

# Filtering and transformation parameters TODO ins smoothing routine
#min_bin_size = snakemake.params.min_bin_size
#min_days_span = snakemake.params.min_days_span
#max_days_span = snakemake.params.max_days_span



# read snv file
print('*'*80 + "\n")
print("Read snv file " + snv_file + "\n")
print('*'*80 + "\n")
seq_info_short_table = io.read_and_extract_snv_file(snv_file)

# read masked positions from vcf file or positions flag if present
if masking_file:
    print('*'*80 + "\n")
    print("Read masked positions " + masking_file+ "\n")
    print('*'*80 + "\n")
   # masked_positions = io.read_vcf_and_extract_masked_positions(masking_file)

masked_positions = []
if masked_positions_str:
    print('*'*80 + "\n")
    print("Parse masked positions\n")
    print('*'*80 + "\n")
    masked_positions = filt.get_positions_from_string(masked_positions_str)

# filtering
print('*'*80 + "\n")
print("Filtering SNV positions\n")
print('*'*80 + "\n")
#seq_info_short_table['snvs_unfiltered'] = seq_info_short_table['snvs']
seq_info_short_table['snvs'] = filt.filter_snvs(seq_info_short_table['snvs'],
                 freq_threshold=freq_cutoff,   
                 masked_positions=masked_positions,
                 remove_indels=False, remove_all_insertions=True)

# binning
print('*'*80 + "\n")
print("***  Binning and phi calculation\n")
print('*'*80 + "\n")
phi_per_bin_table = bin.calculate_phi_per_bin(seq_info_short_table,
                          days_per_bin=days_per_bin,
                          seqs_per_bin=seq_per_bin)

# sequence statistics
print('*'*80 + "\n")
print("*** Infer sequence statistics per day \n")
print('*'*80 + "\n")
seq_info_perDay_table = si.get_seq_info_per_day(seq_info_short_table)

# write tables
print('*'*80 + "\n")
print("***  Write result tables into "+ result_path +"\n")
print('*'*80 + "\n")
io.write_phi_per_bin_table(phi_per_bin_table, path = result_path, suffix=suffix)
io.write_sequence_info_per_day_table(seq_info_perDay_table, path = result_path, suffix=suffix)
