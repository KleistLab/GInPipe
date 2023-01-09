#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 3 2022

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
from ginpipepy.read_masking_vcf import VCFreader

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
        binnings .append(file)
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

# Masking parameters
#masking_file = args.mask_vcf

# Collect all trajectories to merge into one data set
bin_merging_data = []

for folder in binnings:
    print("Current binning folder: %s" % folder)
    # List all folders in a directory
    binnings_dir = "%s/%s" % (bins_dir, folder)

    thetas_comb = []
    variance_comb = []
    variance_size_comb = []
    num_seq_comb = []
    origins_comb = []
    num_mut_comb = []

    list_folders = os.listdir(binnings_dir)
    list_folders.sort()
    list_folders = list_folders[1:len(list_folders)]

    # List header files
    list_headers = [x for x in os.listdir(headers_dir) if not x.startswith('.')]
    headers = []
    for file in list_headers:
        if file.startswith("header"):
            headers.append(file)
    headers.sort()

    for binning_folder in list_folders:
        print("Current position binning folder: %s" % binning_folder)
        # List all folders in a directory
        pos_binnings_dir = "%s/%s/%s" % (bins_dir, folder, binning_folder)
        list_pos_folders = os.listdir(pos_binnings_dir)
        files = list(filter(lambda x: x.endswith("bam"), os.listdir(pos_binnings_dir)))
        #print(list_pos_folders)

        seq_list_base_complete = []
        seq_list_pairs_complete = []

        for filename in files:
            path = files_dir = "%s/%s/%s/%s" % (bins_dir, folder, binning_folder, filename)
            sam_to_fp = SAMtoFP(path,reference,refname)
            seq_lbase, lref, mutants_pairs_list = sam_to_fp.write_fp()
            seq_list_base_complete.append(seq_lbase)
            seq_list_pairs_complete.append(mutants_pairs_list)
        #print(seq_list_base_complete)
        print(" Done.")
        print(" Filtering variants with sequence number cutoff: <=%d" % freqCutoff)
        filt = parameterEstimation(seq_list_base_complete,reference,freqCutoff)
        mut_proportion, filtered_seqset = filt.run()
        print(" Done.")
        #print(len(mut_proportion))
        #print(len(filtered_seqset))

        # Compute optimal metric paramteres for each individual bin
        print(" Computing optimal metric parameters...")
        #Theta from origins - MLE
        analyze = analyzeTrajectory(filtered_seqset, mut_proportion, [1,1,1,1,1], '')
        thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyze_bins_mle()
        print(thetas)
        print(" Done.")
        print(" Writing binning table...")

        # Combine the theta estimates; different options
        # 1. Max
        thetas_max = max(thetas)
        max_ind = thetas.index(thetas_max)
        # Append mean theta estimate over sequence position windows
        # to combine with other temporal estimates
        thetas_comb.append(thetas_max)
        variance_comb.append(variance[max_ind])
        variance_size_comb.append(variance_size[max_ind])
        num_seq_comb.append(num_seqs[max_ind])
        origins_comb.append(origins[max_ind])
        num_mut_comb.append(num_mut[max_ind])

    name_table = "table_%s_phi_estimate_var_from_size.tsv" % (folder)
    table_path = "%s/%s" % ("reconstructed_bins", name_table)
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["date","value","variance_size","variance_mle","num_seqs","times",
                         "num_days_per_bin","haplotypes", "num_mut"])
        for i in range(len(thetas_comb)):
            writer.writerow([thetas_comb[i], variance_size_comb[i], variance_comb[i], num_seq_comb[i],
                            origins_comb[i], num_mut_comb[i]])
    print(" Done.\n")
    # Write to merged bins dataset
    for i, date in enumerate(thetas_comb):
        bin_merging_data.append((thetas_comb[i],
                                variance_size_comb[i],
                                variance_comb[i],
                                num_seq_comb[i],
                                folder,
                                origins_comb[i],
                                num_mut_comb[i]))
print('-'*80)
print(" Making the final results table...")
# Make merged dataset sorted by date
bin_merging_data_ = sorted(bin_merging_data)
times = np.arange(0, len(bin_merging_data_))
name_table = "%s/%s" % ("reconstructed_bins", "table_merged_phi_estimates_var_from_size.tsv")

table_path = "%s" % (name_table)

with open(table_path, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["t","phi","variance_size","variance","num_seq","binning_folder",
                     "origins", "num_mut"])
    for i in range(len(times)):
        # If bin size==1/variance is bigger or equal to min_bin_size
        # if maxsize, minsize satisfied
        #if bin_merging_data_[i][2]<=1/int(min_bin_size) and bin_merging_data_[i][6]>=min_days_span and bin_merging_data_[i][6]<=max_days_span:
        writer.writerow([times[i],
                        bin_merging_data_[i][0],
                        bin_merging_data_[i][1],
                        bin_merging_data_[i][2],
                        bin_merging_data_[i][3],
                        bin_merging_data_[i][4],
                        bin_merging_data_[i][5],
                        bin_merging_data_[i][6]])

print(" Done.")


