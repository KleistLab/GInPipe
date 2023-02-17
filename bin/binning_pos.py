#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import math
from Bio import SeqIO


class BinningPosition:
    def __init__(self, table, qc, ref, dir):
        self.table = table
        self.mean_length = int(qc)
        self.ref = ref
        self.ref_length = 0
        self.ref_seq_list = []
        self.name = ''
        self.dir = dir

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

    def _inferMeanFragmentLength(self):
        # If fastQC output cannot be used
        return 0

    def bin(self):
        """ Write position bin files containing reads from temporal bins
        """
        # Preset bin size and temp files
        # Infer pair length
        self._getReferenceInfo()

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
        
                        
                                

