#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert BAM file into dna profile table and write table.
"""


import bam.bam_to_fingerprints as bam
import utils.io_routines as io

##################################################
### Read parameters
##################################################

bam_file = snakemake.input.bam
ref_file = snakemake.input.ref
result_file =  snakemake.output.snv

print('-'*80 + "\n")
print("Convert BAM to DNA profile table\n")
print('-'*80 + "\n\n")

##################################################
### Convert BAM
##################################################
bam_converter =  bam.SAMtoFP(bam_file,ref_file)

snv_table = bam_converter.get_snv_table()

##################################################
### Write table
##################################################
io.write_table(snv_table, result_file)

