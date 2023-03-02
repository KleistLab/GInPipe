#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import math
import os
from Bio import SeqIO
from matplotlib import pyplot as plt


class BinningPosition:
    def __init__(self, table, qc, ref, dir, bed_file):
        self.table = table
        self.mean_length = int(qc)
        self.ref = ref
        self.ref_length = 0
        self.ref_seq_list = []
        self.name = ''
        self.dir = dir
        self.bed = bed_file

    def _getReferenceInfo(self):
        """Get reference name, sequence and length
        """
        fasta_ref = SeqIO.parse(open(self.ref),'fasta')
        name = ''
        ref_sequence = ''
        for fasta in fasta_ref:
            name, ref_sequence = fasta.id, str(fasta.seq)
        list_sequence = list(ref_sequence.upper())
        self.ref_length = len(list_sequence)
        self.ref_seq_list = list_sequence
        # TODO: name without ">"
        self.name = name

    def _extractBedFile(self):
        t = pd.read_table(self.bed, delimiter="\t", names=['chrom','chromStart','chromEnd','name','primerPool','strand'])
        t_left = t.loc[t['strand'] == '+']
        return t_left

    def _getBinningScheme(self):
        """Record read start positions to select binning windows start
        """
        t = pd.read_table(self.table, delimiter="\t")
        start_positions = t['start_pair']
        start_values, start_counts = np.unique(start_positions, return_counts=True)
        #peaks = find_peaks(start_counts,rel_height=0.2)
        print(start_values)
        print(start_counts)
        #plt.plot(start_values,start_counts)
        #plt.savefig('start_position_distribution.pdf')
        #end_positions = t['end_pair']
        #end_values, end_counts = np.unique(end_positions, return_counts=True)
        return 0

    def _inferMeanFragmentLength(self):
        # If fastQC output cannot be used
        return 0
    
    def _writeFile(self, folder, filename, table):
        outfile = '%s/%s' % (folder, filename)
        table.to_csv(outfile, sep='\t', index=False, mode='w')
    
    '''
    # 
            # Check if current day bigger than bin borders
            bins_bool = [day >= b for b in bins_ranges]
            # Sum  bool list to get bin index
            bin_ind = sum(bins_bool)-1
            bin_index.append(bin_ind)
    '''
    def bin(self):
        """ Write position bin files containing reads from temporal bins
        For primer based binning:
        prim = [prim1, prim2, prim3, ..., primN]
        bool start_pair > prim:
            bool_prim = [True, True, False, ..., False]
        find ind of last True --> index of bin
        """
        # Preset bin size and temp files
        # Infer pair length
        self._getReferenceInfo()
        bed_table = self._extractBedFile()
        # Try binning based on primer starts

        bins_fbed = np.array(bed_table['chromStart'])
        print(bins_fbed)
        temp_bin = pd.read_table(self.table, delimiter="\t")
        bin_index = []

        for index, line in temp_bin.iterrows():
            start = line['start_pair']
            end = line['end_pair']
            # Check if end is within current amplicon area
            bins_bool = [end > p for p in bins_fbed]
            bin_start = sum(bins_bool)-2
            bin_value = bins_fbed[bin_start]
            bin_index.append(bin_value)

        temp_bin["bin_index"] = bin_index

        table_grouped = temp_bin.groupby(['bin_index'])
        keys = table_grouped.groups.keys()
        for i,ampl in enumerate(keys):
            #ind_str = "%s" % i
            bin_table = table_grouped.get_group(ampl)
            filename = "%s_%s.tsv" % ('bin_amplicon',ampl)
            self._writeFile(self.dir ,filename, bin_table)

        return 0

    def bin_st(self):
        """ Write position bin files containing reads from temporal bins
        For primer based binning:
        prim = [prim1, prim2, prim3, ..., primN]
        bool start_pair > prim:
            bool_prim = [True, True, False, ..., False]
        find ind of last True --> index of bin
        """
        # Preset bin size and temp files
        # Infer pair length
        self._getReferenceInfo()
        bed_table = self._extractBedFile()
        # Binning based on primer starts
        bins_fbed = bed_table['chromStart']

        num_bins = math.ceil(self.ref_length/self.mean_length)
        for i in range(num_bins):
            # Make a temporal position bin file to write lines to
            filename = "%s/%s_%s.tsv" % (self.dir,'bin',i)
            with open(self.table, 'r') as infile:
                header = infile.readline()
                next(infile) 
                with open(filename, "w+") as outfile:
                    outfile.write(header)
                    for line in infile:
                        # Split line and get start and stop
                        linesplit = line.split('\t')
                        # name fp1 fp2 start_pair end_pair date
                        # Position in alignment
                        ref_start = int(linesplit[3])
                        ref_end = int(linesplit[4])
                        # 2. By placing read into a bin if its bigger chunk falls into a bin
                        #left_overh = mean_length*i-ref_start # -: inside bin; +: outside bin
                        #right_overh = ref_end-mean_length*(i+1) # -: inside bin; +: ouside bin
                        start_overlap = max(ref_start,self.mean_length*i)
                        end_overlap = min(ref_end,self.mean_length*(i+1)) # If it begins before window -> 
                        len_overlap = end_overlap-start_overlap      # will be smaller than start_overlap
                        query_length = ref_end-ref_start+1
                        if len_overlap > 0:
                            if len_overlap > math.floor(query_length*0.5):
                                # Write read in the corresponding read BAM file
                                outfile.write(line)
        return 0
        
        
                        
                                

