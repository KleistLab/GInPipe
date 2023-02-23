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
import re
import pysam

from bam_to_fingerprints_paired import SAMtoFP
from ginpipepy.modular_theta_from_dict import analyzeTrajectory
from ginpipepy.parameter_est import parameterEstimation


# Get console input
parser = argparse.ArgumentParser()
parser.add_argument("bam", help="BAM directory")
parser.add_argument("ref_fasta", help="Reference FASTA file")
parser.add_argument("region", help="Name of reference region that BAM is aligned to")
#parser.add_argument("mask_vcf", help="VCF file with positions to be masked")
#parser.add_argument("out_dir", help="Output directory")
args = parser.parse_args()

# Output path
bam_file = args.bam

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

# Masking parameters
masking_file = ''

# Getting the mutated positions for each sequence in compressed FP format
print(" Recording mutant positions from CIGAR strings...")
seq_list_pairs_complete = []
path = "%s" % (bam_file)
# Index BAM file
pysam.index(path)
sam_to_fp = SAMtoFP(path,reference,refname)
fp1, fp2, names, start_pair, end_pair = sam_to_fp.write_fp()
print(len(fp1))
print(" Done.")
# Write results
print(" Writing mutant positions from CIGAR strings as fingerprints into file:")
name_table = "table_merged_phi_estimates_var_from_size.tsv"
print("         %s" % name_table)
# Write fingerprint table
table_path = name_table
with open(table_path, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
                            #quotechar='', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["name","fp1","fp2","start_pair","end_pair"])
    for i in range(len(names)):
        writer.writerow([names[i],
                        fp1[i],
                        fp2[i],
                        start_pair[i],
                        end_pair[i]])
print(" Done.")
