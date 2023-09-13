#bam_to_fingerprints.py

"""
Tranlsate CIGAR strings table of dates and SNVs.
"""


import pysam
import random
import numpy as np
from Bio import SeqIO
import pandas as pd
import re
from warnings import warn
import datetime

class SAMtoFP:
    """
    SAMtoFP class.

    Tranlsate CIGAR strings to SNV strings  
    """

    def __init__(self, filename, reffile):
        """
        Turn CIGAR strings in SAM files into sequence fingerprints.

        :param filename: path to SAM/BAM file
        :type filename: str
        :param reffile: path to reference FASTA file
        :type reffile: str
        """
        self.filename = filename
        self.reference = reffile
        self.ref_name = self.get_ref_name(reffile) 

    def get_ref_name(self, reffile):
        with open(str(reffile), "r") as file:
            header = file.readline()
            refname = header.strip(">")
            refname_full = refname.strip("\n")
            split_header = refname_full.split()
            #Minimap2 header
            refname = split_header[0]
            return refname

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

    def _date_from_header(self, name):
        
        date_format = r"[|][\d]{4}-[\d]{2}-[\d]{2}"
        dt = re.search(date_format, name)
        if dt:
            dt = dt[0][1:]
        return dt

    def _cigar_to_fp(self, seq, cigar, start):
        """
        Translate CIGAR strings to sequence fingerprints in format.

            Ref_basePositionMutant_base (e.g. T33G)
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
        #lref = len(ref_seq)
        # Counter reference
        counter = start
        # Counter query
        counter_q = 0
        # Mutants list
        #mutants_base = []
        # Mutants count
        #mutants_pos = 0
        #mutants_pairs_list = []

        snv_list = []

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
                                snv_list.append(ref_seq[counter]+alt_rec +bcall)
                                #mutants_base.append(alt_rec+'>'+bcall)
                                #mutants_pairs_list.append((counter+1,bcall))
                                #mutants_pos += 1
                        else:
                            snv_list.append(ref_seq[counter]+alt_rec +str(alt_base))
                            #mutants_base.append(alt_rec+'>'+str(alt_base))
                            #mutants_pairs_list.append((counter+1,str(alt_base)))
                            #mutants_pos += 1
                    counter += 1
                    counter_q += 1
        # Merge mutant fingerprintts to one string with positions separated by "-"
        snv_str = ''
        if snv_list!=[]:
            snv_str = " ".join(snv_list)

        return snv_str
        #return mutantsstrbase, lref, mutants_pairs_list
    

    def get_snv_table(self):
        """
        TODO describe
        """

        dates= []
        snvs = []
        file = pysam.AlignmentFile(self.filename)

        for read in file.fetch(self.ref_name):
            name = read.query_name
            seq = str(read.query_sequence)
            start_ = read.get_reference_positions()
            start = start_[0]
            # Trim with cigar string
            cigar = read.cigartuples
            date = self._date_from_header(name)
            snv_str = self._cigar_to_fp(seq, cigar, start)
            try:
                datetime.date.fromisoformat(date)
            except ValueError:
                warn("Skipping sequence " + name + ".\nNo date with format yyyy-mm-dd found in header.\n")
            else:
                dates.append(date)
                snvs.append(snv_str)
        
        # sort table by date
        snv_table = pd.DataFrame({'date': dates,
                                'dna_profile': snvs})

        snv_table = snv_table.sort_values(by="date")
        return snv_table

    #TODO weg
    # def write_fp(self):
    #     """
    #     Write sequence fingerprints from CIGAR format in BAM file.

    #     :returns sequences_list_base: fingerprints per bin per sequence as one string
    #     :returns sequence_pos_list: fingerprints per bin per sequence as list of strings
    #     :rtype: list
    #     :returns sequence_pair_list: fingerprints per bin per sequence as pairs of
    #         sequence positions and mutant bases
    #     :rtype: list
    #     :returns lref: length of reference sequence
    #     :rtype: int
    #     """
    #     # Write positions with mutant base as string
    #     sequences_list_base = []
    #     # Write positions with mutant base as pair
    #     sequence_pair_list = []

    #     file = pysam.AlignmentFile(self.filename)

    #     for read in file.fetch(self.ref_name):
    #         name = read.query_name
    #         seq = str(read.query_sequence)
    #         start_ = read.get_reference_positions()
    #         start = start_[0]
    #         # Trim with cigar string
    #         cigar = read.cigartuples

    #         mutants_string, lref, mutants_pairs_list = self._cigar_to_fp(seq, cigar, start, name)

    #         sequences_list_base.append((name,mutants_string))
    #         sequence_pair_list.append(mutants_pairs_list)


    #     return sequences_list_base, lref, sequence_pair_list
