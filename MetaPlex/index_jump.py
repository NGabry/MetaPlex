import pathlib
import math
import sys
import shutil
import tempfile
import os
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins.demux.visualizers import summarize


def calculate(demultiplexed_seqs, sample_map, calibrator_tag_pairs):
    '''
    Calculate Index Jump Rate based off calibrator Tags

    demultiplexed_seqs: path to demultiplexed QIIME2 qza file of data type SampleData[SequencesWithQuality]

    sample_map: path to tab delimited QIIME2 sample map file

    calibrator_tag_pairs: pairs of calibrator tags, input as either a tuple, or list of tuples, with
                          each index being a 2 digit zero padded string.
                          ex: [('01', '11')] or [('01', '11'), ('02','12')]

    Return: Integer value of maximum number of index jumps (false reads) expected in a single sample
    Output: Expected_False_Reads_Per_Index.csv file containing number of false reads expected in EACH sample
            log.txt file containing summary statistics
    '''

    def extract_tsvs(viz, dest):
        '''
        Extracts useful TSV from QIIME artifact files
        '''
        with tempfile.TemporaryDirectory() as temp:
            viz.export_data(temp)
            temp_pathlib = pathlib.Path(temp)
            for file in temp_pathlib.iterdir():
                if file.suffix == '.tsv':
                    shutil.copy(file, dest)

    # Load in sample_map as pandas df
    meta_df = pd.read_csv(sample_map, sep='\t')

    # Load in demultiplexed sequences from path and Summarize
    demux = Artifact.load(demultiplexed_seqs)
    print('Demultiplexed sequences loaded... Extracting per sample summaries to tsvs directory \n . \n . \n .')
    demux_summary = summarize(demux)
    demux_viz = demux_summary.visualization

    # Make directory, extract tsvs, read in per-sample tsv, set as df
    os.makedirs('tsvs', exist_ok=True)
    extract_tsvs(demux_viz, 'tsvs')
    per_sample_df = pd.read_csv("./tsvs/per-sample-fastq-counts.tsv", sep="\t")
    df = per_sample_df

    # Establishing variables and altering dataframes
    # Create FWD index column
    df['Fwd']     = df['sample ID'].str[1:3]

    # Create REV index column
    df['Rev']     = df['sample ID'].str[4:6]
    df            = df.sort_values(['sample ID'])
    df.reset_index(inplace=True)

    # Add in True_False coding
    df['True_False'] = meta_df['True_False']

    if calibrator_tag_pairs is None:
        print('No calibrator tags given. Calculating index jump rate based on 0s in Sample Map.\n '
              '*NOTE* This calculation method is less accurate than when using calibrator tags!\n'
              'Jump rate estimates will be lower than the true rate and may result in retaining\n'
              'false reads within the sample pool.')

        all_false_reads = df[df["True_False"] == 0]
        summed_false = all_false_reads['forward sequence count'].sum()
        jump_rate = (summed_false/(df['forward sequence count'].sum()))*2

    else:
        # Error flagging for Calibrator Tag input: Scripting / Jupyter-Notebook
        print('Validating Calibrator Tag input... \n . \n . \n .')
        if not all(isinstance(item, tuple) for item in calibrator_tag_pairs):
            raise TypeError(
                'Calibrator Tag Index identifiers are not input in proper format. Should be List of Tuples.')

        # Error flagging for Calibrator Tag input: Command Line input
        for i, j in calibrator_tag_pairs:
            if type(i) != str or type(j) != str:
                raise TypeError('Calibrator Tag Index identifiers are not of proper type. Should be type String.')

            if len(i) != 2 or len(j) != 2:
                raise ValueError('Calibrator Tag Index identifiers are not of proper length. Should be 2 digit inputs.')

        # Calculating the jump_rate of the sequencing run based off of indicated calibrator tags
        individual_jump_rates = []
        for fwd, rev in calibrator_tag_pairs:

            fwd_df          = df[df['Fwd'] == fwd]
            rev_df          = df[df['Rev'] == rev]

            # Validity check for calibrator tags
            fwd_True       = fwd_df[fwd_df['True_False'] == 1]
            assert fwd_True['True_False'].sum() == 1, \
                                                      f'Specified tag pair {fwd, rev} does not meet requirements for ' \
                                                      f'usage as calibrator tags. Should only have one true sequence ' \
                                                      f'specified in sample map.'
            rev_True       = rev_df[rev_df['True_False'] == 1]
            assert rev_True['True_False'].sum() == 1, \
                                                      f'Specified tag pair {fwd, rev} does not meet requirements for ' \
                                                      f'usage as calibrator tags. Should only have one true sequence ' \
                                                      f'specified in sample map.'

            print(f'All checks passed! Calculating Index Jump rate based off calibrator tags {fwd, rev}')
            # Calculations
            fwd_False       = fwd_df[fwd_df['True_False'] == 0]
            fwd_jump_rate   = (fwd_False['forward sequence count'].sum())/(fwd_df['forward sequence count'].sum())
            rev_False       = rev_df[rev_df['True_False'] == 0]
            rev_jump_rate   = (rev_False['forward sequence count'].sum())/(rev_df['forward sequence count'].sum())

            individual_jump_rates.append(fwd_jump_rate)
            individual_jump_rates.append(rev_jump_rate)

        # Calculate Average Jump rate for ALL calibrator tags
        jump_series = pd.Series(individual_jump_rates)
        jump_rate = jump_series.sum()/len(jump_series)

    # Preparing lists for calculations
    Indexes       = df['sample ID'].tolist()
    Seq_counts    = df['forward sequence count'].tolist()

    Fwds          = df['Fwd'].tolist()
    Fwds          = list(set(Fwds))
    Fwds.sort()

    Revs          = df['Rev'].tolist()
    Revs          = list(set(Revs))
    Revs.sort()

    # Calculating each Indexes expected # of False reads
    indexes_exp_false = []
    total_read_count = df['forward sequence count'].sum()

    count = 0
    for i in Fwds:
        temp_df = df[df['Fwd'] == i]

        # No. of reads with i FWD index
        i_seqs = temp_df['forward sequence count'].sum()
        i_seqs_no_j = i_seqs - Seq_counts[count]

        for j in Revs:
            temp_df2 = df[df['Rev'] == j]
            # No. of reads with j REV index
            j_seqs = temp_df2['forward sequence count'].sum()
            j_seqs_no_i = j_seqs-Seq_counts[count]

            # No. of j reads with an expected false i index
            i_jumps = i_seqs * jump_rate
            i_to_j = i_jumps * (j_seqs_no_i / total_read_count)

            # No. of i reads with an expected false j index
            j_jumps = j_seqs * jump_rate
            j_to_i = j_jumps * (i_seqs_no_j / total_read_count)

            # Sum for total i-j pair false read count
            indiv_tag_jump_count = math.ceil((i_to_j + j_to_i))

            # Append to list for data table output
            to_add = (Indexes[count], indiv_tag_jump_count)
            indexes_exp_false.append(to_add)
            count += 1

    # Export csv with expected false reads per sample
    out_df = pd.DataFrame(indexes_exp_false, columns=['SampleIndex', 'Expected False Reads'])
    print('Exporting recommended per-sample filtering levels to Expected_False_Reads_Per_Index.csv')
    out_df.to_csv('Expected_False_Reads_Per_Index.csv', index=False)

    # Text output to log.txt of some summary statistics
    max_IJR = out_df['Expected False Reads'].max()
    false_read_count = out_df['Expected False Reads'].sum()

    with open('log.txt', 'w') as out:
        print(f'Calculated Average jump rate: {jump_rate}', file=out)
        print(f'Total number of false reads: {false_read_count} / {total_read_count}', file=out)
        print(f'Total percent of false reads: {(false_read_count/total_read_count)*100:.3f}%', file=out)
        print(f'Maximum number of false reads expected in a single sample: {max_IJR}', file=out)

    return max_IJR


def main():
    if len(sys.argv) < 4:
        calculate(sys.argv[1], sys.argv[2], calibrator_tag_pairs=None)

    else:
        c_tags_list = []
        for i in range(3, len(sys.argv)):
            c_tag_pair = (tuple(sys.argv[i].split(',')))
            c_tags_list.append(c_tag_pair)
            calculate(sys.argv[1], sys.argv[2], c_tags_list)


if __name__ == '__main__':
    main()
