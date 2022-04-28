import os
import subprocess
import pandas as pd
import sys


def remultiplex(sequenceFile, indexFile):
    '''
    Takes dual-indexed reads, trims the 5' and 3' ends of the reads past the indexes, and moves the 3' index to
    immediately follow the 5' index (i.e. ['MultiplexedSingleEndBarcodeInSequence'] format)

    sequenceFile : path to raw sequence file of type .fastq, .fastq.gz, or .bam

    indexFile    : path to .csv containing all the index tag sequences that are present in the sequencing pool. This .csv
                   should be formatted as specified below (See indexes.csv for reference)

    Return : None

    Output : remultiplexed_seqs.fastq.gz
    '''
    if sequenceFile.endswith('.bam'):
        convert = str('samtools fastq '+sequenceFile+' > raw_seqs.fastq')
        subprocess.run(convert, shell=True)
        sequenceFile = './raw_seqs.fastq'

    else:
        assert sequenceFile.endswith('.fastq') or sequenceFile.endswith('.gz'), 'Sequence file not of proper format. Should be .bam or .fastq'

    # Read in indexes csv, and convert Forward indexes to dictionary
    df = pd.read_csv(indexFile, dtype=str)
    df['ID_ori'] = (df['orientation'] + df['ID'])
    f_df = df[df['orientation'] == 'F']
    fwd_indexes = dict(zip(f_df.ID_ori, f_df.seq))

    # Organizing by Forward Index
    f = open('fwd_sort.txt', 'w')
    for k in fwd_indexes:
        fwd_sort = str('cutadapt -g '+fwd_indexes[k]+' --discard-untrimmed -e 0 -j 0 --overlap 10 -o '
                   +k+'.fastq '+sequenceFile)
        subprocess.run(fwd_sort, shell=True, stdout=f)

    os.system('find . -type f -size 0 -delete')

    # Convert Reverse barcodes to dictionary
    r_df = df[df['orientation'] == 'R']
    rev_indexes = dict(zip(r_df.ID_ori, r_df.seq))

    fastqs = []
    for file in os.listdir('.'):
        if len(file) == 9 and file.endswith('fastq'):
            fastqs.append(file)

    # Iterating over each Forward Index to then organize by Reverse index
    f = open('rev_sort.txt', 'w')
    for fastq in fastqs:
        for k in rev_indexes:
            rev_sort = str('cutadapt -a '+rev_indexes[k]+'$ --discard-untrimmed -e 0 -j 0 --overlap 10 -o '
                           +fastq[0:3]+k+'.fastq '+fastq)
            subprocess.run(rev_sort, shell=True, stdout=f)

    os.system('find . -type f -size 0 -delete')

    # Create dictionary for merged Fwd and Rev idnexes
    merged_indexes = {}
    for k_a, v_a in fwd_indexes.items():
        for k_b, v_b in rev_indexes.items():
            merged_indexes[k_a + k_b] = v_a + v_b

    fastqs = []
    for fastq in os.listdir('.'):
        if len(fastq) == 12 and fastq.endswith('fastq'):
            fastqs.append(fastq)

    # Prepending each read with the appropriate merged index as indicated by file name
    fastqs.sort()
    for fastq in fastqs:
        for k, v in merged_indexes.items():
            if fastq[0:6] == k:
                with open(fastq, 'r') as f, open('remultiplexed_seqs.fastq', 'a+') as out:
                    count = 0
                    for idx, line in enumerate(f.read().splitlines()):
                        count += 1
                        if count == 2:
                            print(v+line, file=out)
                        elif count == 4:
                            print(''.join(['I']*len(v))+line, file=out)
                        else:
                            print(line, file=out)
                        if count == 4:
                            count = 0
            else:
                continue

    os.system('gzip remultiplexed_seqs.fastq')

    os.system('mkdir fastqs')
    os.system('mv *.fastq fastqs| mv fwd_sort.txt fastqs| mv rev_sort.txt fastqs')


def main():
    remultiplex(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
