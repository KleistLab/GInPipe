#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 12:16:30 2020

@authors: Yannick Duport, Maria Trofimova
"""

"""
Reads in fasta-file and tsv-file, containing metadata, and combines the two by adding the date from the meta-data to
the name and id of the fasta-entries
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import sys
import argparse


parser=argparse.ArgumentParser(
    description='''Merge date from metadata table with FASTA headers in GISAID format. ''')
parser.add_argument('--i', type=str, default=None, help='Input FASTA file')
parser.add_argument('--m', type=str, default=None, help='Input meta file')
parser.add_argument('--mc', type=str, default=None, help='Column name by which metadata will be merged with FASTA headers - column containing full headers from the FASTA file')
parser.add_argument('--dc', type=str, default=None, help='Column name containing date that will be added to FASTA headers')
parser.add_argument('--o', type=str, default=None, help='Output FASTA file')
args=parser.parse_args()

fasta_file = vars(args)['i']
meta_file = vars(args)['m']
merge_column = vars(args)['mc']
date_column = vars(args)['dc']
out_file = vars(args)['o']

# read csv and return dict with key=strain, value=strain|date
meta_dict = {}
with open (meta_file, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        meta_dict[row[merge_column]] = f"{row[merge_column]}|{row[date_column]}"

# remove XX from dates (e.g 2020-03-XX -> 2020-03)
for strain, date in meta_dict.items():
    if date.endswith("-XX"):
        meta_dict[strain] = date[:-3]

# read in fasta-file, and create new seq-records
new_records = []
seq_records = SeqIO.parse(open(fasta_file), 'fasta')
for record in seq_records:
    name = meta_dict[record.name]
    new_records.append(SeqRecord(record.seq, id=name, name=name, description=''))

# write new seq-records into new fasta file
SeqIO.write(new_records, out_file, 'fasta')

