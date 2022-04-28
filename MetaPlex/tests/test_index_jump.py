import os
import subprocess
import unittest
from MetaPlex import index_jump


class IndexJumpTests(unittest.TestCase):

    def test_single_cal(self):
        self.val = index_jump.calculate(demultiplexed_seqs='data/demux_seqs/demultiplexed_seqs.qza',
                                        sample_map='data/Sample_Map.txt',
                                        calibrator_tag_pairs=[('01', '11')])

        self.assertEqual(self.val, 5)

        # Clean up
        os.system('rm Expected_False_Reads_Per_Index.csv | rm -r tsvs | rm log.txt')

    def test_no_cal(self):
        self.val = index_jump.calculate(demultiplexed_seqs='data/demux_seqs/demultiplexed_seqs.qza',
                                        sample_map='data/Sample_Map.txt',
                                        calibrator_tag_pairs=None)

        self.assertEqual(self.val, 1)
        # Clean up
        os.system('rm Expected_False_Reads_Per_Index.csv | rm -r tsvs | rm log.txt')


if __name__ == '__main__':
    unittest.main()
