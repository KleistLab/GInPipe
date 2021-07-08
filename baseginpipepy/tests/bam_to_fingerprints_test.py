import unittest
from ginpipepy.bam_to_fingerprints import SAMtoFP

"""
Created on Tue Jul 6 2021

@authors: Maria Trofimova
"""

class testFP(test_bam,test_reference,test_reference_name):

    def __init__(self, bam_path, reference)
        self.bam_path = bam_path
        self.reference = reference

    def testFingerprint(self):
        #Reference name
        with open(str(reference), "r") as file:
            header = file.readline()
        refname = header.strip(">")
        refname = refname.strip("\n")

        sam_to_fp = SAMtoFP(self.path,self.reference,refname)
        seq_lbase, seq_posbase, lref, mutants_pairs_list = sam_to_fp.writeFP()
        fp_result = '45>A-51>A-56>A-85>G-124>A-126>C-194>A'
        print("Correct sequence fingerpint:")
        print(fp_result)
        print("Sequence fingerprint in test sequence")
        print(seq_lbase[0].value)
