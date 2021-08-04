#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:24:55 2020

@author: Yannick Duport, Maria Trofimova
"""


from pathlib import Path
from ginpipepy.sam_to_bins_modular import SAM
import math
import pandas as pd
import os
import csv
import subprocess

file_name = snakemake.input.bam[0]
stats_name = snakemake.input.stats[0]
RESULT_PATH = Path(os.getcwd()) / "results"
bins_dir = RESULT_PATH / "bins"
meta_dir = RESULT_PATH / "meta"
SAM_PATH = file_name # Main big SAM or BAM file
WORKING_PATH = RESULT_PATH

num_per_bin = snakemake.params.eq_num
# If no explicit binning desired - use default fractions
if not num_per_bin:
    # Read table containing file statistics created by samtools
    stats = pd.read_table(str(stats_name),delimiter='\t',header=None)
    num_reads = stats.iloc[0,2]
    num_per_bin.append(int(math.floor(num_reads*0.02)))
    num_per_bin.append(int(math.floor(num_reads*0.05)))
    num_per_bin.append(int(math.floor(num_reads*0.07)))
    print("No sizes for binning by equal size provided, binning will use the\n following default size fractions:\n")
    print("      * 2%%, or %i sequences" % int(num_per_bin[0]))
    print("      * 5%%, or %i sequences" % int(num_per_bin[1]))
    print("      * 7%%, or %i sequences" % int(num_per_bin[2]))
    print('-'*80)

days_per_bin = snakemake.params.eq_days
param1 = num_per_bin[0]
param2 = days_per_bin[0]

reffile_dir = str(snakemake.params.reference)
with open(str(reffile_dir), "r") as file:
    header = file.readline()
param3 = header.strip(">")
param3header = param3.strip("\n")
#Split header on whitespaces -- if whitespaces present the forst string will be taken
split_header = param3header.split()
#Minimap2 header
mm_header = split_header[0]
print('Reference sequence header full: %s' % param3header)
print('Reference sequence header used: %s' % mm_header)
sam = SAM(SAM_PATH, bins_dir, meta_dir, param1, param2, mm_header)

print('-' * 80)
print("Bins of equal Size")
print('-' * 80)
for param in num_per_bin:
    sam.eq_num_param = param
    sam.bin_eq_size_names()

for param in days_per_bin:
    sam.eq_days_param = param
    print('\n')
    print('-'*80)
    print("Bins with equal number of days")
    print('-'*80)
    sam.bin_eq_days()
    print('\n')
    print('-'*80)
    #print("Bins with equal number of days (fuzzy)")
    #print('-'*80)
    #sam.bin_eq_days(fuzzy=True)

print('\n')
print('-'*80)
print("Bins by calendar weeks")
print('-'*80)
sam.bin_cal_week()

'''
Create tsv files containing all bin file names per binning folder
'''
# List binnings
bins = os.listdir(str(bins_dir))
name_file = "list_of_binnings.tsv"
file_path = "%s/%s" % (bins_dir, name_file)
with open(file_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(bins)):
            writer.writerow([bins[i]])
