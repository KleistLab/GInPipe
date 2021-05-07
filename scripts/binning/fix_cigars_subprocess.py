#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:15:49 2020

@author: Maria Trofimova
"""

from pathlib import Path
import os
import csv
import subprocess

FILEPATH = Path(__file__).parent
list_file_name = snakemake.input[0]
tool_path = snakemake.params[0]
ref = snakemake.params[2]
FILEPATH_interm = FILEPATH.parent
ref_path = FILEPATH_interm.parent / "consensus" / ref
ref_path_str = str(ref_path)+'.fasta'
results_dir = FILEPATH_interm.parent / "results"
WORKING_PATH = FILEPATH_interm.parent

list_folders = os.listdir(str(results_dir)+"/bins")
list_folders_filter = list(filter(lambda x: not x.startswith('.') and not x.endswith('.tsv'), list_folders))
print(list_folders_filter)

for folder in list_folders_filter:
    print(folder)
    list_files = os.listdir(str(results_dir)+'/bins/'+folder)
    list_files_filter = list(filter(lambda x: x.endswith('.bam'), list_files))
    list_files_filter.sort()
    if not os.path.exists(str(results_dir)+'/fixed_cigars_bins/'+folder):
        os.makedirs(str(results_dir)+'/fixed_cigars_bins/'+folder)
    for file in list_files_filter:
        infile_path = str(results_dir)+'/bins/'+folder+'/'+file
        outsam_path = str(results_dir)+'/fixed_cigars_bins/'+folder+'/'+file[:-4]+'.sam'
        outfixsam_path = str(results_dir)+'/fixed_cigars_bins/'+folder+'/'+file[:-4]+'_fixcigar.sam'
        outfixbam_path = str(results_dir)+'/fixed_cigars_bins/'+folder+'/'+file[:-4]+'_fixcigar.bam'
        subprocess.call("samtools view -h -o %s %s" % (outsam_path,infile_path), shell=True)
        subprocess.call(	"java -jar %s -r %s %s > %s" % (tool_path,ref_path_str,outsam_path,outfixsam_path), shell=True)
        subprocess.call("rm %s" % outsam_path, shell=True)
        subprocess.call("samtools view -S -b %s > %s" % (outfixsam_path,outfixbam_path), shell=True)
        subprocess.call("rm %s" % outfixsam_path, shell=True)
        subprocess.call("samtools index %s" % (outfixbam_path), shell=True)

    name_file = "list_of_files.tsv"
    file_path = str(results_dir)+'/fixed_cigars_bins/'+folder+"/"+name_file
    with open(file_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(list_files_filter)):
            writer.writerow([list_files_filter[i]])


name_file = "list_of_binnings.tsv"
file_path = str(results_dir)+'/fixed_cigars_bins/'+name_file
with open(file_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(list_folders_filter)):
            writer.writerow([list_folders_filter[i]])
