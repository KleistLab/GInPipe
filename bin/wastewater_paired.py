#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 14:23:20 2022

@author: mariatrofimova
"""

import pysam
import random
import itertools
import numpy as np

import math
import pandas as pd
import os

from Bio import SeqIO
from pathlib import Path

class wastewaterReadsToSeq:
    def __init__(self, bam_file, region, ref, out_dir, seed):
        self.bam = pysam.AlignmentFile(bam_file, "rb")
        self.region = region
        self.ref = ref
        self.fastqc = ''
        self.out_dir = out_dir
        # For testing
        self.seed = seed
        self.ref_length = 0
        self.dates = []
        self.ref_seq_list = []
        self.read_counts = []

    def _writePseudoSequence(self, sequence, ind, date):
        """Write pseudoreads to a fasta file
        """
         # Write read to fasta file
        out = self.out_dir
        if not os.path.exists(out):
            os.makedirs(out)
        filename = 'ww_simulated.fasta'
        # Make sequence name
        # TODO: make GISAID-like format for date in sequence headers
        with open(os.path.join(out, filename), 'a') as f:
            f.write('>WW_sim_sequence_'+str(ind)+'|'+str(date))
            f.write("\n")
            # Write sequence
            f.write(sequence)
            f.write("\n")
            f.close()
        return filename

    def _createTempBinFile(self, ind):
        filename = "bin_file_"+str(ind)+".tmp"
        filedir = Path(self.out_dir)
        filepath = filedir / filename
        os.makedirs(filedir, exist_ok=True)
        return filepath

    def _createBinFile(self, ind):
        filename = "bin_file_"+str(ind)+".bam"
        filedir = Path(self.out_dir) 
        filepath = filedir / filename
        os.makedirs(filedir, exist_ok=True)
        return filepath

    def _getReferenceInfo(self):
        """Get reference sequence name, sequence and length
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

    def _inferMeanLength(self):
        """Infer length without FastQC stats
        """
        reads_iter = self.bam.fetch(self.region)
        read_lengths = []
        for r in reads_iter:
            l = r.tlen()
            read_lengths.append(l)
        return int(math.ceil(np.mean(read_lengths))*1.5), sum(read_lengths)

    def _getFastqcStats(self):
        """ Summarize the FastQC stats 
        """
        # Prior to that:
            #  awk '/>>Sequence Length Distribution/{flag=1; next} />>END_MODULE/{flag=0} flag' file
        rfile = pd.read_csv(self.fastqc,sep='\t')
        # Mean read length
        count = rfile['Count']
        # Total number of bases estimate
        length = rfile['Length']
        sizes = [float(str(l).split('-')[1]) for l in length]
        lmean = math.ceil(np.mean(sizes))
        # Upper estimate of total bases count
        total_bases = sum(sizes*count)
        return lmean, total_bases

    def _writeBinFiles(self):
        """ Write bin files containing reads from bins
        """
        # Bin size from FastQC stats
        mean_length, total_bases = self._inferMeanLength()
        # TODO: Get info about the reference sequence
        num_bins = math.ceil(self.ref_length/mean_length)
        # Count reads in bins
        read_counts = []
        for i in range(num_bins-1):
            print("    Writing bin #",i)
            # Fetch the iterator over reads in BAM file
            reads_iter = self.bam.fetch(self.region)
            filepath = self._createTempBinFile(i)
            # Read count
            count = 0
            with pysam.AlignmentFile(str(filepath), "w", header=self.bam.header) as temp_bin_file:
                for read in reads_iter:
                    # Read tid to zero; see issue https://github.com/pysam-developers/pysam/issues/509
                    read.tid = 0
                    # Header with date
                    header = read.query_name
                    self.dates.append(header.split("|")[-1])
                    # Position in alignment
                    ref_start = read.reference_start
                    # Create the bin file via pysam
                    # 1. By checking if start position falls into the bin
                    if (ref_start in range(mean_length*i,mean_length*(i+1))):
                        # Write read in the corresponding read BAM file
                        temp_bin_file.write(read)
                        count += 1
                        if read.next_reference_start >= ref_start:
                            next                    
            # Add count to list
            read_counts.append(count)
            # SAM to BAM
            infile = pysam.AlignmentFile(filepath, "r")
            bin_file = self._createBinFile(i)
            with  pysam.AlignmentFile(bin_file, "wb", template=infile) as bamfile:
                for s in infile:
                    bamfile.write(s)
            os.remove(filepath)
            # Index the bin BAM file 
            pysam.index(str(bin_file))
        # TODO: nicer message printing, maybe via some functions?
        # Add read counts 
        self.read_counts = read_counts
        print("   Done writing bins.")
        return num_bins

    def _writeBinFiles2(self):
        """ Write bin files containing reads from bins
        """
        # Bin size from FastQC stats
        mean_length, total_bases = self._inferMeanLength()
        # TODO: Get info about the reference sequence
        num_bins = math.ceil(self.ref_length/mean_length)
        # Count reads in bins
        read_counts = []
        for i in range(num_bins-1):
            print("    Writing bin #",i)
            # Fetch the iterator over reads in BAM file
            reads_iter = self.bam.fetch(self.region)
            filepath = self._createTempBinFile(i)
            # Read count
            count = 0
            with pysam.AlignmentFile(str(filepath), "w", header=self.bam.header) as temp_bin_file:
                for read in reads_iter:
                    # Read tid to zero; see issue https://github.com/pysam-developers/pysam/issues/509
                    read.tid = 0
                    # Header with date
                    header = read.query_name
                    self.dates.append(header.split("|")[-1])
                    # Position in alignment
                    ref_start = read.reference_start
                    ref_end = read.reference_start+read.query_length
                    # Create the bin file via pysam
                    # 2. By placing read into a bin if its bigger chunk falls into a bin
                    #left_overh = mean_length*i-ref_start # -: inside bin; +: outside bin
                    #right_overh = ref_end-mean_length*(i+1) # -: inside bin; +: ouside bin
                    start_overlap = max(ref_start,mean_length*i)
                    end_overlap = min(ref_end,mean_length*(i+1)) # If it begins before window -> 
                    len_overlap = end_overlap-start_overlap      # will be smaller than start_overlap
                    if len_overlap > 0:
                        if len_overlap > read.query_length or len_overlap >= mean_length:
                            # Write read in the corresponding read BAM file
                            temp_bin_file.write(read)
                            count += 1
                            if read.next_reference_start >= ref_start:
                                next     
                    #if ((mean_length*(i+1)-ref_start)/read.query_length > 0.5):
                    #    temp_bin_file.write(read)
                    #    count += 1
                    #    if read.next_reference_start >= ref_start:
                    #        next
            # Add count to list
            read_counts.append(count)
            # SAM to BAM
            infile = pysam.AlignmentFile(filepath, "r")
            bin_file = self._createBinFile(i)
            with  pysam.AlignmentFile(bin_file, "wb", template=infile) as bamfile:
                for s in infile:
                    bamfile.write(s)
            os.remove(filepath)
            # Index the bin BAM file 
            pysam.index(str(bin_file))
        # TODO: nicer message printing, maybe via some functions?
        # Add read counts 
        self.read_counts = read_counts
        print("   Done writing bins.")
        return num_bins

    def _sampleSequence(self):
        """ Sample one full-length sequence
        """
        # List the bins directory
        bins_dir = Path(self.out_dir)
        filenames = [_ for _ in os.listdir(str(bins_dir)) if _.endswith(".bam")]
        # Record sampled mutations
        mutations_dict = {}
        # Combine length of reads sampled
        length = 0

        for i in range(len(filenames)):
            # TODO check if fetch region works on truncated files
            # Filepath 
            filepath = str(bins_dir)+'/'+filenames[i]
            # Index the bin BAM file 
            pysam.index(str(filepath))
            # Get reads iterator 
            bin_file = pysam.AlignmentFile(filepath, "rb")
            reads_iter = bin_file.fetch(self.region)
            for r in reads_iter:
                # Generate a random float to sample
                rand = random.uniform(0,1)
                if self.read_counts[i] > 0:
                    if rand < 1/self.read_counts[i]:
                        # Get the read length
                        l = r.infer_query_length()
                        length = length+l
                        # Get aligning positions and reference base
                        aligned_pairs = r.get_aligned_pairs(with_seq=True)
                        # Get query sequence
                        query = list(r.query_sequence)
                        # Add mutations to dictionary
                        for i in range(len(aligned_pairs)):
                            # If there is a mismatch
                            if (aligned_pairs[i][0]!=None and aligned_pairs[i][1]!=None):
                                if (aligned_pairs[i][2]!=query[aligned_pairs[i][0]]) and (aligned_pairs[i][2]!=None):
                                    # Record position in frequency array
                                    # Tuple: (reference_pos,reference_base,query_base)
                                    mutations_dict[aligned_pairs[i][1]] = query[aligned_pairs[i][0]].upper()
                        break
        return mutations_dict, length

    def _sampleSequence(self):
        """ Sample one full-length sequence
        """
        # List the bins directory
        bins_dir = Path(self.out_dir)
        filenames = [_ for _ in os.listdir(str(bins_dir)) if _.endswith(".bam")]
        # Record sampled mutations
        mutations_dict = {}
        # Combine length of reads sampled
        length = 0

        for i in range(len(filenames)):
            # TODO check if fetch region works on truncated files
            # Filepath 
            filepath = str(bins_dir)+'/'+filenames[i]
            # Index the bin BAM file 
            pysam.index(str(filepath))
            # Get reads iterator 
            bin_file = pysam.AlignmentFile(filepath, "rb")
            reads_iter = bin_file.fetch(self.region)
            for r in reads_iter:
                # Generate a random float to sample
                rand = random.uniform(0,1)
                if self.read_counts[i] > 0:
                    if rand < 1/self.read_counts[i]:
                        # Get the read length
                        l = r.infer_query_length()
                        length = length+l
                        # Get aligning positions and reference base
                        aligned_pairs = r.get_aligned_pairs(with_seq=True)
                        # Get query sequence
                        query = list(r.query_sequence)
                        # Add mutations to dictionary
                        for i in range(len(aligned_pairs)):
                            # If there is a mismatch
                            if (aligned_pairs[i][0]!=None and aligned_pairs[i][1]!=None):
                                if (aligned_pairs[i][2]!=query[aligned_pairs[i][0]]) and (aligned_pairs[i][2]!=None):
                                    # Record position in frequency array
                                    # Tuple: (reference_pos,reference_base,query_base)
                                    mutations_dict[aligned_pairs[i][1]] = query[aligned_pairs[i][0]].upper()
                        break
        return mutations_dict, length

    def main(self):
        """ Main 
        """
        # Random seed 
        if (self.seed!=0):
            random.seed(self.seed)

        # Record reference information
        self._getReferenceInfo()
        # Maximum length of all reads sampled
        mean_length, total_bases = self._inferMeanLength()
        # Mutation dictionaries, one per simulated sequence
        sequence_dicts = []
        # Total length of simulated sequences - stop if exceeds total_bases
        curr_length = 0
        # Write bins
        num_bins = self._writeBinFiles2()

    def mainFull(self):
        """ Main 
        """
        # Random seed 
        if (self.seed!=0):
            random.seed(self.seed)

        # Record reference information
        self._getReferenceInfo()
        # Maximum length of all reads sampled
        mean_length, total_bases = self._inferMeanLength()
        # Mutation dictionaries, one per simulated sequence
        sequence_dicts = []
        # Total length of simulated sequences - stop if exceeds total_bases
        curr_length = 0
        # Write bins
        num_bins = self._writeBinFiles()
        print("Sampling sequences...")
        ind = 0
        print("Total bases: ",total_bases)
        while curr_length < total_bases:
            mutations_dict, length = self._sampleSequence()
            sequence_dicts.append(mutations_dict)
            curr_length = curr_length+length
            ind += 1
            if (ind % 1000)==0:
                print("     Sampled",ind,"reads")
                print("            and",100*round(curr_length/total_bases,3),"% of bases in sample")
        # Get the reference sequence
        print("Writing FASTA file")
        list_sequence = self.ref_seq_list
        # Make new sequences by adding mutations to reference 
        fin_fasta = ''
        for i,rd in enumerate(sequence_dicts):
            sim_sequence_list = list_sequence
            for pos, base in rd.items():
                sim_sequence_list[pos] = base
            sim_sequence = ''.join(sim_sequence_list)
            datedf = pd.DataFrame(data={"date":self.dates})
            date_ = str(datedf["date"].astype('datetime64[ns]').quantile(0.5, interpolation="midpoint"))
            date = date_.split(" ")[0]
            fin_fasta = self._writePseudoSequence(sim_sequence,i,date)
        print("Done")
        #return fin_fasta





    