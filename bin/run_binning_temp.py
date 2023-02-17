#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from binning_temp import BinningTemporal


parser = argparse.ArgumentParser()
parser.add_argument("table", help="TSV file")
#parser.add_argument("mode", choices=['size', 'days'], help="Binning by size or by number of days")
#parser.add_argument("size", help="Parameter for binning mode")
#parser.add_argument("work_dir", help="Working directory")

args = parser.parse_args()
os.mkdir('bins')

""" temp_binning = BinningTemporal(args.table, args.size)
if args.mode == "size":
    temp_run = temp_binning.binBySize()
elif args.mode == "days":
    temp_run = temp_binning.binByDays()
else:
    raise argparse.ArgumentTypeError('Invalid mode parameter. Valid options are size and days.')
 """

# Binning by days only
sizes = [1]

for size in sizes:
    temp_binning = BinningTemporal(args.table, size)
    temp_run = temp_binning.binByDays()
