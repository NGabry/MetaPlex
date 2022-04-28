import os
import subprocess
import unittest
from MetaPlex import remultiplexing


class RemuxTest(unittest.TestCase):
    def test_valid(self):
        remultiplexing.remultiplex(sequenceFile='data/raw_fastq/raw_seqs.fastq.gz',
                                   indexFile='data/indexes.csv')
        out = subprocess.Popen(["wc", "-l", "remultiplexed_seqs.fastq.gz"],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)

        stdout, stderr = out.communicate()
        self.result = int(stdout.split()[0].decode("utf-8"))
        self.assertEqual(self.result, 304686)

        # Clean up
        os.system('rm remultiplexed_seqs.fastq.gz')
        os.system('rm -r fastqs')


if __name__ == '__main__':
    unittest.main()
