#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 4 2022

@author: Maria Trofimova
"""

import argparse
import os
from pathlib import Path
import numpy as np
import math
from math import gamma
import csv
import pandas as pd
import datetime
import random
from collections import Counter

from ginpipepy.bam_to_fingerprints import SAMtoFP
from ginpipepy.modular_theta_from_dict import analyzeTrajectory
from ginpipepy.parameter_est import parameterEstimation
from ginpipepy.read_masking_vcf import VCFreader

# Get console input
parser = argparse.ArgumentParser()
parser.add_argument("bins", help="Bins directory")
parser.add_argument("ref_fasta", help="Reference FASTA file")
parser.add_argument("region", help="Name of reference region that BAM is aligned to")
parser.add_argument("group", help="Grouping string: type of binning")
parser.add_argument("freq_cutoff", help="Mutant frequency cutoff")
parser.add_argument("mask_vcf", help="VCF file with positions to be masked")
parser.add_argument("out", help="Output directory") # Might not need this
#parser.add_argument("out_dir", help="Output directory")
args = parser.parse_args()


# Output path
out_dir = args.out
bins_dir = args.bins
#header_dir = args.bins
list_binnings = os.listdir(str(bins_dir))
binnings = []
for file in list_binnings:
    if not file.startswith(".") and not file.endswith(".tsv"):
        binnings.append(file)
binnings.sort()

# Reference location
reference = str(args.ref_fasta)
region = args.region

# Grouping variable - suffix for all files
file_suffix = args.group

# Base frequency cutoff
freqCutoff = args.freq_cutoff

# Masking parameters
masking_file = args.mask_vcf

for folder in binnings:
    print("Current binning folder: %s" % folder)
    # List all folders in a directory
    binnings_dir = "%s/%s" % (bins_dir, folder)
    list_folders = os.listdir(binnings_dir)
    list_folders.sort()

    #
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

        # Count partial haplotypes in filtered_seqset
        part_haps = Counter(filtered_seqset)
        C = 100
        N = len(filtered_seqset)
        x = 1
        # Likelihood function
        L1 = gamma(N+1)/math.prod(gamma(x+1), start=1)
        L2 = gamma(sum(C*part_haps.keys))/gamma(N+math.sum(C*part_haps.keys))
        L3 = math.prod(gamma(x+C*part_haps.keys)/gamma(C*part_haps.keys))
        L = L1*L2*L3
        



