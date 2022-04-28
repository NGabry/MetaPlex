import sys
import biom
import pandas as pd
from qiime2 import Artifact


def per_sample_filter(feature_table, filtering_integer):
    '''
    Filters reads out of a QIIME2 feature table according to a minimum read count requirement *per sample*

    feature_table    : path to QIIME2 feature table

    filtering_integer: Either an integer for even filtering across samples, or path to the
                       Expected_False_Reads_Per_Index.csv output by index_jump.py

    Return: Frequency filtered QIIME2 feature table of type FeatureTable[Frequency]

    Output: 'freq_filt_table.qza' QIIME2 artifact of type FeatureTable[Frequency]
    '''

    # Importing from a biom table
    try:
        working_table = biom.load_table(feature_table)
    except:
        print('Could not import from biom table, trying from qza...')
        try:
            # Importing from imported qza
            qza_table = Artifact.load(feature_table)
            working_table = qza_table.view(biom.Table)
        except:
            print('Could not import from qiime artifact, check to ensure input file is appropriate format')
            exit()

    # Convert table to df
    biom_df = working_table.to_dataframe()

    # Filtering if integer is specified
    if type(filtering_integer) == int:

        print('File loaded! Filtering at '+str(filtering_integer)+' occurrences per feature')

        # Filter out any feature from a sample if it has less than X reads,
        biom_df = biom_df._get_numeric_data()
        biom_dense = biom_df.sparse.to_dense()

        biom_dense[biom_dense < filtering_integer] = 0

        # Translate the filtered table back into a qza
        biom_dense = biom_dense.transpose()
        freq_filt_table = Artifact.import_data("FeatureTable[Frequency]", biom_dense)

        # Export the new table as a qza for further use in Qiime2
        freq_filt_table.save('freq_filt_table.qza')

        return freq_filt_table

    # Filtering based on per-sample CSV
    else:
        print(f'File loaded! Filtering according to per-sample values specified in {filtering_integer}')
        per_sample_table = pd.read_csv(filtering_integer, index_col='SampleIndex')
        per_sample_table_trans = per_sample_table.transpose()

        biom_df = biom_df._get_numeric_data()
        biom_dense = biom_df.sparse.to_dense()

        for k, v in per_sample_table_trans.iteritems():
            for column, values in biom_dense.iteritems():
                if k == column:
                    working_values = values.tolist()
                    edited_values = [0 if x <= int(v) else x for x in working_values]
                    biom_dense[column] = edited_values

        # Translate the filtered table back into a qza
        biom_dense = biom_dense.transpose()
        freq_filt_table = Artifact.import_data("FeatureTable[Frequency]", biom_dense)

        # Export the new table as a qza for further use in Qiime2
        freq_filt_table.save('freq_filt_table.qza')

        return freq_filt_table


def main():
    if sys.argv[2].endswith('.csv'):
        per_sample_filter(sys.argv[1], sys.argv[2])
    else:
        per_sample_filter(sys.argv[1], int(sys.argv[2]))


if __name__ == '__main__':
    main()
