#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:24:55 2020

@author: Yannick Duport, Maria Trofimova
"""


from pathlib import Path
from sam_to_bins_modular import SAM
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
param3 = param3.strip("\n")
sam = SAM(SAM_PATH, bins_dir, meta_dir, param1, param2, param3)

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
Create tsv files containing all bin file names per binning folder and erase empty bins
'''
bins = os.listdir(str(bins_dir))
bins_filter = list(filter(lambda x: not x.startswith('.'), bins))

for folder in bins_filter:
    files = os.listdir(str(bins_dir)+"/"+folder)
    headers = list(filter(lambda x: x.startswith('header'), files))
    files_filter = list(filter(lambda x: x.startswith('bin'), files))
    ranges_filter = list(filter(lambda x: x.startswith('range'), files))
    files_filter.sort()
    headers.sort()
    ranges_filter.sort()
    for i,hfile in enumerate(headers):
        table = pd.read_table(str(bins_dir)+"/"+folder+"/"+hfile,delimiter='\t',header=0)
        if table.empty:
            bam_name = str(bins_dir)+"/"+folder+"/"+files_filter[i]
            header_name = str(bins_dir)+"/"+folder+"/"+headers[i]
            range_name = str(bins_dir)+"/"+folder+"/"+ranges_filter[i]
            subprocess.call("rm %s" % bam_name, shell=True)
            subprocess.call("rm %s" % header_name, shell=True)
            subprocess.call("rm %s" % range_name, shell=True)
        else:
            bam_name = str(bins_dir)+"/"+folder+"/"+files_filter[i]
            subprocess.call("samtools index %s" % bam_name, shell=True)
    files = os.listdir(str(bins_dir)+"/"+folder)
    files_filter = list(filter(lambda x: x.endswith('.bam'), files))
    files_filter.sort()
    name_file = "list_of_files.tsv"
    file_path = str(bins_dir)+"/"+folder+"/"+name_file
    with open(file_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(files_filter)):
            writer.writerow([files_filter[i]])

name_file = "list_of_binnings.tsv"
file_path = str(bins_dir)+'/'+name_file
with open(file_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(bins_filter)):
            writer.writerow([bins_filter[i]])
