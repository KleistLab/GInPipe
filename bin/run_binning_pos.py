#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from binning_pos import BinningPosition


parser = argparse.ArgumentParser()
parser.add_argument("dir", help="Directory of temporal binnings")
#parser.add_argument("size", help="Mean read length")
parser.add_argument("ref", help="Reference file")
#parser.add_argument("work_dir", help="Working directory")

args = parser.parse_args()

# List binning directory
bin_modes = os.listdir(args.dir)
os.mkdir('pos_bins')

for bin_mode in bin_modes:
    pos_dir = '%s/%s' % ('pos_bins', str(bin_mode))
    dir = '%s/%s' % (args.dir, str(bin_mode))
    print(pos_dir)
    print(dir)
    os.mkdir(pos_dir)

    for table in os.listdir(dir):   
        table_dir = '%s/%s' % (dir, table)
        print(table_dir)
        out_table_dir = '%s/%s' % (pos_dir, table[0:-4])
        print(out_table_dir)
        os.mkdir(out_table_dir)
        pos_binning = BinningPosition(table_dir, 385, args.ref, out_table_dir)
        pos_binning.bin()
