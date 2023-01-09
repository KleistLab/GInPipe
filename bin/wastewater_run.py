#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:51:17 2022

@author: mariatrofimova
"""



import argparse
import os
import pandas as pd
from pathlib import Path
from datetime import datetime
from wastewater import wastewaterReadsToSeq


# Get console input
parser = argparse.ArgumentParser()
parser.add_argument("bins", help="Bins directory")
parser.add_argument("ref_fasta", help="Reference FASTA file")
parser.add_argument("region", help="Name of reference region that BAM is aligned to")
parser.add_argument("num_repeats",help="Number of times the daaset should be run (for simulation testing)")
#parser.add_argument("out_dir", help="Output directory")
args = parser.parse_args()


# Run with input parameters
bins_path = args.bins
bins_dir = Path(bins_path)
print(bins_dir)
bins_dir_list__ = os.listdir(str(bins_dir))
bins_dir_list = list(filter(lambda x: not x.startswith(('bin','list')), bins_dir_list__)) 
#bins_dir_list = list(filter(lambda x: not x.startswith('list'), bins_dir_list_))

num_repeats = int(args.num_repeats)


def main(bins_dir_list,
         ref_fasta,
         region):
    # 
    for folder in bins_dir_list:
        bin_path = bins_dir / folder
        bins = os.listdir(str(bin_path))
        bins_filter_ = list(filter(lambda x: not x.endswith('bai'), bins))
        bins_filter = sorted(list(filter(lambda x: not x.endswith('tsv'), bins_filter_)))
        
        out_dir__ = Path(os.getcwd()) / 'sorted_binned'
        out_dir_ = out_dir__ / folder
        for i, bin_ in enumerate(bins_filter):
            bam = bin_path / bin_
            # Split folder per one temporal bin file
            out_dir = out_dir_ / bin_.rsplit('.', 1)[0]
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            ww = wastewaterReadsToSeq(bam_file=str(bam),
                    ref=ref_fasta,
                    region=region,
                    out_dir=str(out_dir),
                    seed=0)
            ww.main()
    return 0
    
result = main(bins_dir_list,
         args.ref_fasta,
         args.region)