#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 12:38:18 2021

@author: mariatrofimova
"""

from pathlib import Path
import pandas as pd
import os
import csv
import subprocess

'''
Create tsv file containing all bin file names per binning folder without empty bins
'''
bins_path = snakemake.input[0]
bins_ppath = Path(bins_path)
bins_dir = bins_ppath.parent.absolute()
print(bins_dir)
bins = os.listdir(str(bins_dir))
print(bins)
bins_filter = list(filter(lambda x: not x.startswith('list'), bins))

for folder in bins_filter:
    files = os.listdir("%s/%s" % (bins_dir, folder))
    headers = list(filter(lambda x: x.startswith('header'), files))
    files_filter = list(filter(lambda x: x.startswith('bin'), files))
    ranges_filter = list(filter(lambda x: x.startswith('range'), files))
    files_filter.sort()
    headers.sort()
    ranges_filter.sort()
    for i,hfile in enumerate(headers):
        table = pd.read_table("%s/%s/%s" % (bins_dir, folder, hfile),delimiter='\t',header=0)
        if table.empty:
            bam_name = "%s/%s/%s" % (bins_dir, folder, files_filter[i])
            header_name = "%s/%s/%s" % (bins_dir, folder, headers[i])
            range_name = "%s/%s/%s" % (bins_dir, folder, ranges_filter[i])
            os.remove(bam_name)
            os.remove(header_name)
            os.remove(range_name)
        else:
            bam_name = "%s/%s/%s" % (bins_dir, folder, files_filter[i])
            subprocess.check_call("samtools index %s" % bam_name, shell=True)

name_file = "list_of_binnings_filtered.tsv"
file_path = "%s/%s" % (bins_dir, name_file)
with open(file_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(bins_filter)):
            writer.writerow([bins_filter[i]])
