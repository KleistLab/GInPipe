#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:05:41 2020

@author: Maria Trofimova
"""
import pysam
import random
import numpy as np
from Bio import SeqIO


class SAMtoFP:

    def __init__(self, filename, reffile, refname):
        self.filename = filename
        self.reference = reffile
        self.ref_name = refname
        
    def _ambiguousDict(self, base, refbase):
        """
        Calls ambiguous bases if reference base is not included in the possible calls
        :param base: ambiguous base
        :param refbase: reference base
        :return variant: what base to call
        """
        ambdict = {'W':['T','A'],
                   'S':['C','G'],
                   'M':['A','C'],
                   'K':['G','T'],
                   'R':['A','G'],
                   'Y':['C','T'],
                   'B':['C','G','T'],
                   'D':['A','G','T'],
                   'H':['A','C','T'],
                   'V':['A','G','C'],
                   'N':['A','C','G','T'],
                   'Z':[]}
        # Get bases in corresponding ambiguous base list
        possvar = ambdict[base]
        variant = base
        # If reference base is in list - do not call a variant
        if refbase in possvar:
            variant = refbase
        # Else - randomly decide which base to call based on ambiguous list
        else:
            u = random.uniform(0,1)
            prob_base = np.full(len(possvar),1/len(possvar))
            cprob_base = np.hstack((0,np.cumsum(prob_base)))
            for i in range(1,len(cprob_base)):
                if cprob_base[i-1]<u<cprob_base[i]:
                    variant = possvar[i-1]

        return variant
            
    
    def _cigarToFP(self, seq, cigar, start, name):
        """
        Translates CIGAR strings to sequence fingerprints in format 
        Position>Mutant_base
        CIGAR string signatures are taken in form of Pysam cigartuple in format:
            (operation,length)
        Key list:
            M (mis-match) = 0
            I (insertion) = 1
            D (deletion) = 2
            N (reference skip) = 3
            S (soft clip) = 4
            H (hard clip) = 5
            P (pad) = 6
            = (match) = 7
            X (mismatch) = 8
            B (back) = 9
        :param seq: query sequence
        :param cigar: cigar string as pysam cigar table
        :param start: position on query where the query starts
        :param name: name of query sequence
        :return mutantsstrbase: fingerprints per sequence as one string
        :return mutants_pos_list: fingerprints per sequence as list of strings
        :return mutants_pairs_list: fingerprints per sequence as pairs of 
            sequence positions and mutant bases
        """
        # SeqIO the reference file
        ref_fasta = SeqIO.parse(open(self.reference),'fasta')
        ref_seq = ''
        for fasta in ref_fasta:
            ref_seq = str(fasta.seq)
        lref = len(ref_seq)
        # Counter reference
        counter = start
        # Counter query
        counter_q = 0
        # Mutants list
        mutants_base = []
        # Mutants count
        mutants_pos = 0
        mutants_pos_list = []
        mutants_pairs_list = []
       
        # Ambiguous bases list
        amblist = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N', 'Z']
        for i in range(len(cigar)):
            if cigar[i][0]==0:
                counter += cigar[i][1]
                #if counter > 21544:
                #    break
                counter_q += cigar[i][1]
            elif cigar[i][0]==1:
                counter_q += cigar[i][1]
            elif cigar[i][0]==2:
                counter += cigar[i][1]
                #if counter > 21544:
                #   break
            elif cigar[i][0]==3:
                counter += cigar[i][1]
                #if counter > 21544:
                 #   break
            elif cigar[i][0]==4:
                counter_q += cigar[i][1]
            elif  cigar[i][0]==7:
                counter += cigar[i][1]
                #if counter > 21544:
                #    break
                counter_q += cigar[i][1]
            # Record base that is a mismatch (X) on both query and reference counters
            elif cigar[i][0]==8:
                # Also write the mutation event in the string/fingerprint +ref_seq[counter]+'>'
                for j in range(cigar[i][1]):
                    alt_rec = str(counter+1)
                    alt_base = seq[counter_q]
                    if str(alt_base)!='N':
                        # Record mutant base
                        if str(alt_base) in amblist:
                            # Ambiguous call
                            bcall = self._ambiguousDict(str(alt_base), ref_seq[counter])
                            if bcall==ref_seq[counter]:
                                pass
                            else:
                                mutants_base.append(alt_rec+'>'+bcall)
                                mutants_pos_list.append(counter+1)
                                mutants_pairs_list.append((counter+1,bcall))
                                mutants_pos += 1
                        else:
                            mutants_base.append(alt_rec+'>'+str(alt_base))
                            mutants_pos_list.append(counter+1)
                            mutants_pairs_list.append((counter+1,str(alt_base)))
                            mutants_pos += 1
                    counter += 1
                    counter_q += 1
        # Merge mutant fingerprintts to one string with positions separated by "-"
        mutantsstrbase = ''
        if mutants_base!=[]:
            mutantsstrbase = "-".join(mutants_base)
        
        return mutantsstrbase, mutants_pos_list, lref, mutants_pairs_list


    def writeFP(self):
        """
        Write sequence fingerprints from CIGAR format in BAM file
        :return sequences_list_base: fingerprints per bin per sequence as one string
        :return sequence_pos_list: fingerprints per bin per sequence as list of strings
        :return sequence_pair_list: fingerprints per bin per sequence as pairs of 
            sequence positions and mutant bases
        :return lref: length of reference sequence
        """
        # Write pairs
        sequences_list_base = []
        # Write pure positions
        sequence_pos_list = []
        # Write positions with mutant base
        sequence_pair_list = []
       
        file = pysam.AlignmentFile(self.filename)

        for read in file.fetch(self.ref_name):
            name = read.query_name
            seq = str(read.query_sequence)
            start_ = read.get_reference_positions()
            start = start_[0]
            # Trim with cigar string
            cigar = read.cigartuples
            
            mutants_string, mutants_pos_list, lref, mutants_pairs_list = self._cigarToFP(seq, cigar, start, name)

            sequences_list_base.append((name,mutants_string))
            sequence_pos_list.append((name,mutants_pos_list))
            sequence_pair_list.append(mutants_pairs_list)
            

        return sequences_list_base, sequence_pos_list, lref, sequence_pair_list
