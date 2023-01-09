#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 16:44:06 2020

@author: Maria Trofimova
"""

import argparse
import os
from pathlib import Path
import numpy as np
import math
import csv
import pandas as pd
import datetime
import random


from ginpipepy.bam_to_fingerprints import SAMtoFP
from ginpipepy.modular_theta_from_dict import analyzeTrajectory
from ginpipepy.parameter_est import parameterEstimation


# Get console input
parser = argparse.ArgumentParser()
parser.add_argument("bins", help="Bins directory")
parser.add_argument("header", help="Headers directory")
parser.add_argument("ref_fasta", help="Reference FASTA file")
parser.add_argument("region", help="Name of reference region that BAM is aligned to")
parser.add_argument("freq_cutoff", help="Mutant frequency cutoff")
#parser.add_argument("mask_vcf", help="VCF file with positions tto be masked")
#parser.add_argument("out_dir", help="Output directory")
args = parser.parse_args()

# Placeholders default values - change later
min_bin_size = 1
min_days_span = 1

# Output path
bins_dir = args.bins
headers_dir = args.header
list_binnings = os.listdir(str(bins_dir))
binnings = []
for file in list_binnings:
    if not file.startswith(".") and not file.endswith(".tsv"):
        binnings.append(file)
binnings.sort()


# Reference location
reference = str(args.ref_fasta)

# Name of reference sequence
with open(str(reference), "r") as file:
    header = file.readline()
refname = header.strip(">")
refname_full = refname.strip("\n")
split_header = refname_full.split()
#Minimap2 header
refname = split_header[0]
region = args.region

# Grouping variable - suffix for all files
#file_suffix = args.group

# Base frequency cutoff
freqCutoff = int(args.freq_cutoff)

# Filtering and transformation parameters
min_bin_size = 1
min_days_span = 1
max_days_span = 5

# Masking parameters
masking_file = ''

# Collect all trajectories to merge into one data set
bin_merging_data = []

# List header files
list_headers = [x for x in os.listdir(headers_dir) if not x.startswith('.')]
headers = []
for file in list_headers:
    if file.startswith("header"):
        headers.append(file)
headers.sort()


#List individual binning directories
for folder in binnings:
    print("Current binning folder: %s" % folder)
    # List all files in a directory
    binnings_dir = "%s/%s" % (bins_dir, folder)
    headers_dir = "%s/%s" % (headers_dir, folder)
    list_files = os.listdir(binnings_dir)
    directory = os.listdir(binnings_dir)
    # List BAM files
    files = []
    for file in directory:
        if file.endswith(".bam"):
            files.append(file)
    files.sort()
    # List header files
    list_headers = os.listdir(headers_dir)
    headers = []
    for file in list_headers:
        if file.startswith("header"):
            headers.append(file)
    headers.sort()

    # Initialize result
    # Positions with mutant base as string
    seq_list_base_complete = []
    # Positions with mutant base as pair
    seq_list_pairs_complete = []

    # Getting the mutated positions for each sequence in compressed FP format
    print(" Recording mutant positions from CIGAR strings...")
    for filename in files:
        path = "%s/%s" % (binnings_dir, filename)
        sam_to_fp = SAMtoFP(path,reference,refname)
        seq_lbase, lref, mutants_pairs_list = sam_to_fp.write_fp()
        seq_list_base_complete.append(seq_lbase)
        seq_list_pairs_complete.append(mutants_pairs_list)
    print(" Done.")
    print(" Filtering variants with sequence number cutoff: <=%d" % freqCutoff)
    filt = parameterEstimation(seq_list_base_complete,reference,freqCutoff)
    mut_proportion, filtered_seqset = filt.run()
    print(" Done.")

    # Masking sites specified in optional VCF file
    #if masking_file:
    #    print(" Masking bases based on provided Variant Calling File...")
    #    mask = VCFreader(masking_file, reference, filtered_seqset)
    #    filtered_seqset = mask.mask_bases_in_fp()
    #    print(" Done.")

    # Get the time span and mean date of each bin
    mean_header_bin = []
    # Internal times vector from 0 to N
    times = []
    # Number of days per bin - check compliance with filter
    num_days_per_bin = []
    i = 0
    for header_file in headers:
        # Read header table
        path = "%s/%s" % (headers_dir,header_file)
        table = pd.read_table(path, header=0)
        if not table.empty:
            # List dates column
            dates = table['date'].tolist()
            dates_dt = np.array(dates, dtype='datetime64[s]')
            # Days difference
            delta_days = (max(dates_dt) - min(dates_dt)) / np.timedelta64(1, 'D') + 1
            num_days_per_bin.append(delta_days)
            # Mean date of bin
            mean = (np.array(dates, dtype='datetime64[s]').view('i8').mean().astype('datetime64[s]'))
            mean_header_bin.append(str(mean)[:10])
            times.append(i)
            i += 1

    # Compute opttimal metric paramteres for each individual bin
    print(" Computing optimal metric parameters...")
    #Theta from origins - MLE
    analyze = analyzeTrajectory(filtered_seqset, mut_proportion, num_days_per_bin, '')
    thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyze_bins_mle()
    print(" Done.")
    print(" Writing binning table...")

    # Write bin file
    name_table = "table_%s_%s_phi_estimate_var_from_size.tsv" % ('example',folder)
    table_path = "%s" % (name_table)
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["date","value","variance_size","variance_mle","num_seqs","times",
                         "num_days_per_bin","haplotypes", "num_mut"])
        for i in range(len(thetas)):
            writer.writerow([mean_header_bin[i], thetas[i], variance_size[i], variance[i], num_seqs[i], times[i],
                             num_days_per_bin[i], origins[i], num_mut[i]])
    print(" Done.\n")
    # Write to merged bins dataset
    for i, date in enumerate(mean_header_bin):
        if not i+1 == len(mean_header_bin):
            bin_merging_data.append((mean_header_bin[i],
                                    thetas[i],
                                    variance_size[i],
                                    variance[i],
                                    num_seqs[i],
                                    times[i],
                                    num_days_per_bin[i],
                                    folder,
                                    origins[i],
                                    num_mut[i]))

print('-'*80)
print(" Making the final results table...")
# Make merged dataset sorted by date
bin_merging_data_ = sorted(bin_merging_data)
times = np.arange(0, len(bin_merging_data_))
name_table = "table_merged_phi_estimates_var_from_size.tsv"

table_path = name_table

with open(table_path, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["t","value","variance","meanBinDate","sampleSize",
                     "daysPerBin", "binning", "haplotypes", "numMut"])
    for i in range(len(times)):
        # If bin size==1/variance is bigger or equal to min_bin_size
        # if maxsize, minsize satisfied
        if bin_merging_data_[i][2]<=1/int(min_bin_size) and bin_merging_data_[i][6]>=min_days_span and bin_merging_data_[i][6]<=max_days_span:
            writer.writerow([times[i],
                            bin_merging_data_[i][1],
                            bin_merging_data_[i][2],
                            bin_merging_data_[i][0],
                            bin_merging_data_[i][4],
                            bin_merging_data_[i][6],
                            bin_merging_data_[i][7],
                            bin_merging_data_[i][8],
                            bin_merging_data_[i][9]])
print(" Done.")
