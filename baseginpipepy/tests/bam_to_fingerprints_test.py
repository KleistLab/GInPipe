import unittest
from ginpipepy.bam_to_fingerprints import SAMtoFP

# [('WS|500||7d63c560-6ece-4afe-bc53-80bd827cddd0|2020-01-16', '45>A-51>A-56>A-85>G-124>A-126>C-194>A')]

class TestFP(test_bam,test_reference,test_reference_name):

    def testBinning(self):
        sam_to_fp = SAMtoFP(path,reference,refname)
        seq_lbase, seq_posbase, lref, mutants_pairs_list = sam_to_fp.writeFP()
        fp_result = '45>A-51>A-56>A-85>G-124>A-126>C-194>A'
        self.assertEqual(seq_lbase[0],value, fp_result, "Wrong string")

if __name__ == '__main__':
    unittest.main()
