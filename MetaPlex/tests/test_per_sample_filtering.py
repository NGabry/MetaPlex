import os
import biom
import unittest


from MetaPlex import per_sample_filtering


class PerSampleFilterTest(unittest.TestCase):
    def test_csv(self):
        freq_filt_table = per_sample_filtering.per_sample_filter(feature_table='data/per_sample_artifacts/feature_table.qza',
                                               filtering_integer='data/per_sample_artifacts/Expected_False_Reads_Per_Index.csv')

        calc_table = freq_filt_table.view(biom.Table)
        calc_df = calc_table.to_dataframe()
        self.result = calc_df.to_numpy().sum()
        self.assertEqual(self.result, 369915)

        # Clean up
        os.system('rm freq_filt_table.qza')

    def test_single_int(self):
        freq_filt_table = per_sample_filtering.per_sample_filter(feature_table='data/per_sample_artifacts/feature_table.qza',
                                               filtering_integer=500)

        calc_table = freq_filt_table.view(biom.Table)
        calc_df = calc_table.to_dataframe()
        self.result = calc_df.to_numpy().sum()
        self.assertEqual(self.result, 178357)

        # Clean up
        os.system('rm freq_filt_table.qza')


if __name__ == '__main__':
    unittest.main()
