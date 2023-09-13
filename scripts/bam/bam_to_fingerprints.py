#bam_to_fingerprints.py

"""
Tranlsate CIGAR strings to fingerprint strings.

Created on Mon Jul 13 15:05:41 2020

@author: Maria Trofimova
"""


import pysam
import random
import numpy as np
from Bio import SeqIO

class SAMtoFP:
    """
    SAMtoFP class.

    Tranlsate CIGAR strings to fingerprint strings  
    """

    def __init__(self, filename, reffile, refname):
        """
        Turn CIGAR strings in SAM files into sequence fingerprints.

        :param filename: path to SAM/BAM file
        :type filename: str
        :param reffile: path to reference FASTA file
        :type reffile: str
        :param refname: name of reference sequence
        :type refname: str
        """
        self.filename = filename
        self.reference = reffile
        self.ref_name = refname

    def _ambiguous_dict(self, base, refbase):
        """
        Call ambiguous bases if reference base is not included in the possible calls.

        :param base: ambiguous base
        :type base: str
        :param refbase: reference base
        :type refbase: str
        :returns variant: what base to call
        :rtype: string char
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


    def _cigar_to_fp(self, seq, cigar, start, name):
        """
        Translate CIGAR strings to sequence fingerprints in format.

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
        :type seq: str
        :param cigar: cigar string as pysam cigar table
        :type cigar: list
        :param start: position on query where the query starts
        :type start: int
        :param name: name of query sequence
        :type name: str
        :returns mutantsstrbase: fingerprints per sequence as one string
        :rtype: list
        :returns mutants_pos_list: fingerprints per sequence as list of strings
        :rtype: list
        :returns mutants_pairs_list: fingerprints per sequence as pairs of
            sequence positions and mutant bases
        :rtype: list
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
        mutants_pairs_list = []

        # Ambiguous bases list
        amblist = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N', 'Z']
        for i,cigtuple in enumerate(cigar):
            operation, length = cigtuple
            if operation==0:
                counter += length
                #if counter > 21544:
                #    break
                counter_q += length
            elif operation==1:
                counter_q += length
            elif operation==2:
                counter += length
                #if counter > 21544:
                #   break
            elif operation==3:
                counter += length
                #if counter > 21544:
                 #   break
            elif operation==4:
                counter_q += length
            elif operation==7:
                counter += length
                #if counter > 21544:
                #    break
                counter_q += length
            # Record base that is a mismatch (X) on both query and reference counters
            elif operation==8:
                # Also write the mutation event in the string/fingerprint +ref_seq[counter]+'>'
                for j in range(length):
                    alt_rec = str(counter+1)
                    alt_base = seq[counter_q]
                    if str(alt_base)!='N':
                        # Record mutant base
                        if str(alt_base) in amblist:
                            # Ambiguous call
                            bcall = self._ambiguous_dict(str(alt_base), ref_seq[counter])
                            if bcall==ref_seq[counter]:
                                pass
                            else:
                                mutants_base.append(alt_rec+'>'+bcall)
                                mutants_pairs_list.append((counter+1,bcall))
                                mutants_pos += 1
                        else:
                            mutants_base.append(alt_rec+'>'+str(alt_base))
                            mutants_pairs_list.append((counter+1,str(alt_base)))
                            mutants_pos += 1
                    counter += 1
                    counter_q += 1
        # Merge mutant fingerprintts to one string with positions separated by "-"
        mutantsstrbase = ''
        if mutants_base!=[]:
            mutantsstrbase = "-".join(mutants_base)

        return mutantsstrbase, lref, mutants_pairs_list


    def write_fp(self):
        """
        Write sequence fingerprints from CIGAR format in BAM file.

        :returns sequences_list_base: fingerprints per bin per sequence as one string
        :returns sequence_pos_list: fingerprints per bin per sequence as list of strings
        :rtype: list
        :returns sequence_pair_list: fingerprints per bin per sequence as pairs of
            sequence positions and mutant bases
        :rtype: list
        :returns lref: length of reference sequence
        :rtype: int
        """
        # Write positions with mutant base as string
        sequences_list_base = []
        # Write positions with mutant base as pair
        sequence_pair_list = []

        file = pysam.AlignmentFile(self.filename)

        for read in file.fetch(self.ref_name):
            name = read.query_name
            seq = str(read.query_sequence)
            start_ = read.get_reference_positions()
            start = start_[0]
            # Trim with cigar string
            cigar = read.cigartuples

            mutants_string, lref, mutants_pairs_list = self._cigar_to_fp(seq, cigar, start, name)

            sequences_list_base.append((name,mutants_string))
            sequence_pair_list.append(mutants_pairs_list)


        return sequences_list_base, lref, sequence_pair_list
