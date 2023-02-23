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
import sys
import argparse
#sys.path.append(".") 

parser = argparse.ArgumentParser()
parser.add_argument("bam", help="BAM file")
parser.add_argument("stats", help="Samtools stats file")
parser.add_argument("ref", help="Reference FASTA file")
#parser.add_argument("work_dir", help="Working directory")


args = parser.parse_args()

file_name = args.bam
stats_name = args.stats
RESULT_PATH = Path(os.getcwd())
bins_dir = RESULT_PATH / "bins"
meta_dir = RESULT_PATH / "meta"
SAM_PATH = file_name # Main big SAM or BAM file
WORKING_PATH = os.getcwd()

num_per_bin = []

def main(file_name,
         stats_name,
         RESULT_PATH,
         bins_dir,
         meta_dir,
         SAM_PATH,
         WORKING_PATH,
         num_per_bin):
    
    
    days_per_bin = [1]
    param1 = 100
    param2 = days_per_bin[0]
    
    reffile_dir = str(args.ref)
    with open(str(reffile_dir), "r") as file:
        header = file.readline()
    param3 = header.strip(">")
    param3header = param3.strip("\n")
    #Split header on whitespaces -- if whitespaces present the first string will be taken
    split_header = param3header.split()
    #Minimap2 header
    mm_header = split_header[0]
    print('Reference sequence header full: %s' % param3header)
    print('Reference sequence header used: %s' % mm_header)
    sam = SAM(SAM_PATH, bins_dir, meta_dir, param1, param2, mm_header)
    
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
    
    
    '''
    Create tsv files containing all bin file names per binning folder
    '''
    # List binnings
    # bins = os.listdir(str(bins_dir))
    # name_file = "list_of_binnings.tsv"
    # file_path = "%s/%s" % (RESULT_PATH, name_file)
    # with open(file_path, 'w+', newline='') as csvfile:
    #         writer = csv.writer(csvfile, delimiter='\t',
    #                                 quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #         for i in range(len(bins)):
    #             writer.writerow([bins[i]])
    
    return 0


result = main(file_name,
              stats_name,
              RESULT_PATH,
              bins_dir,
              meta_dir,
              SAM_PATH,
              WORKING_PATH,
              num_per_bin)