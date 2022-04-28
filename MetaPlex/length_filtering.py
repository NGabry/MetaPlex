import os
import sys
import pathlib
import shutil
import tempfile

import pandas as pd
from Bio import SeqIO
from qiime2 import Artifact, Metadata
from qiime2.plugins.feature_table.methods import filter_features, filter_seqs


# Length filter of rep_seqs
def length_filter(feature_table, representative_sequences, length_to_filter):
    '''
    Length filter of QIIME2 feature table and representative sequences

    feature_table: path to QIIME2 qza file of data type FeatureTable[Frequency]

    representative_sequences: path to QIIME2 qza file of data type FeatureData[Sequence]

    length_to_filter: An integer threshold for sequence length. All sequences shorter than specified
                      length will be removed.

    Return: Frequency filtered QIIME2 feature table of type FeatureTable[Frequency]
            Frequency filtered QIIME2 feature table of type FeatureData[Sequence]

    Output: 'length_filt_table.qza' QIIME2 artifact of type FeatureTable[Frequency]
            'length_filt_seqs.qza' QIIME2 artifact of type FeatureData[Sequence]
    '''

    # Pull fasta from rep_seqs to use as metadata
    def extract_fasta(file, dest):
        with tempfile.TemporaryDirectory() as temp:
            file.export_data(temp)
            temp_pathlib = pathlib.Path(temp)
            for item in temp_pathlib.iterdir():
                if item.suffix == '.fasta':
                    shutil.copy(item, dest)

    # QIIME Artifact Loading and error checking
    try:
        qza_table = Artifact.load(feature_table)
    except:
        print('Could not import feature table, check to ensure input file is appropriate format')
        exit()
    try:
        rep_seqs = Artifact.load(representative_sequences)
    except:
        print('Could not import sequences, check to ensure input file is appropriate format')
        exit()

    print(f'QIIME Artifacts loaded! Filtering reads at a length of {length_to_filter}bp')

    os.makedirs('fastas', exist_ok=True)
    extract_fasta(rep_seqs, 'fastas')

    # Parse fasta to get IDs and lengths, then send to df and filter IDs by length
    with open('fastas/dna-sequences.fasta') as fasta_file:
        featureids = []
        lengths = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            featureids.append(seq_record.id)
            lengths.append(len(seq_record.seq))

    rep_seqs_meta = pd.DataFrame({"FeatureID": featureids, "Length": lengths})

    features_to_exclude = rep_seqs_meta[rep_seqs_meta['Length'] < length_to_filter]

    features_to_exclude['FeatureID'].to_csv('Features-to-exclude.csv', index=False)

    exclude = Metadata.load("Features-to-exclude.csv")

    shutil.rmtree('fastas')

    # Filter rep-seqs based on seqs_to_exclude
    rep_seqs_filt = filter_seqs(rep_seqs,
                               metadata = exclude,
                               exclude_ids = True)

    # Filter table based on seqs_to_exclude
    table_filt = filter_features(qza_table,
                                metadata = exclude,
                                exclude_ids = True)

    rep_seqs_filt.filtered_data.save(f'length_filt_seqs_{length_to_filter}.qza')
    table_filt.filtered_table.save(f'length_filt_table_{length_to_filter}.qza')

    return table_filt, rep_seqs_filt


def main():
    length_filter(sys.argv[1], sys.argv[2], int(sys.argv[3]))


if __name__ == '__main__':
    main()

