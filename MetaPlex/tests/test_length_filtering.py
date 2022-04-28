import os
import subprocess
import pathlib
import tempfile
import shutil
import biom
import unittest
from MetaPlex import length_filtering


class LengthFilterTest(unittest.TestCase):
    def test_seqs_filter(self):
        table, seqs = length_filtering.length_filter(feature_table='data/length_artifacts/feature_table.qza',
                                       representative_sequences='data/length_artifacts/rep_seqs.qza',
                                       length_to_filter=150)

        def extract_fasta(file, dest):
            with tempfile.TemporaryDirectory() as temp:
                file.export_data(temp)
                temp_pathlib = pathlib.Path(temp)
                for item in temp_pathlib.iterdir():
                    if item.suffix == '.fasta':
                        shutil.copy(item, dest)

        os.makedirs('fastas', exist_ok=True)
        extract_fasta(seqs.filtered_data, 'fastas')
        out = subprocess.Popen(["wc", "-l", "fastas/dna-sequences.fasta"],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)

        stdout, stderr = out.communicate()
        self.result = int(stdout.split()[0].decode("utf-8"))
        self.assertEqual(self.result, 3662)

        # Clean up
        shutil.rmtree('fastas')

    def test_table_filter(self):
        table, seqs = length_filtering.length_filter(feature_table='data/length_artifacts/feature_table.qza',
                                       representative_sequences='data/length_artifacts/rep_seqs.qza',
                                       length_to_filter=150)

        calc_table = table.filtered_table.view(biom.Table)
        calc_df = calc_table.to_dataframe()
        self.result = calc_df.to_numpy().sum()
        self.assertEqual(self.result, 368114)

        # Clean up
        os.system('rm length_filt_table_150.qza length_filt_seqs_150.qza Features-to-exclude.csv')


if __name__ == '__main__':
    unittest.main()
