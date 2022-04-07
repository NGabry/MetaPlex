# MetaPlex

***MetaPlex*** is a library preparation workflow and read processing toolkit for efficient and accurate COI
metabarcoding on Ion Torrent sequencers. At its core, MetaPlex utilizes multiple pairs of uniquely indexed fusion
primers which contain Ion sequencing adapters, Ion Xpress Barcodes, and the COI ANML primers, allowing for single-step
PCR library preparation of dual-indexed reads.

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/full_read.png?raw=true)

In order to easily process MetaPlex reads in popular analysis platform [QIIME2](https://qiime2.org), MetaPlex provides a
toolkit capable of the following tasks:

* **Remultiplexing** : reorganize dual-indexes to follow ['MultiplexedSingleEndBarcodeInSequence'] format

* **Index jump calculating** : calculate the rate at which index jumps occur per sequencing run, and estimate the number
  of false reads within the total pool as well as each individual sample within the pool

* **Per-sample frequency based filtering** : remove false reads from each sample based off the calculated expectancies
  either provided by the Index Jump calculator, or at user specified levels.

## Installation

All MetaPlex tools function properly when installed in
a [QIIME2 conda environment](https://docs.qiime2.org/2022.2/install/native/#install-qiime-2-within-a-conda-environment) (>
=2021.11)

      conda install -c conda-forge metaplex

# Remultiplexing

Function: takes dual-indexed reads, trims the 5' and 3' ends of the reads past the indexes, and moves the 3' index to
immediately follow the 5' index (i.e. ['MultiplexedSingleEndBarcodeInSequence'] format)

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/standard_dual_index.png?raw=true)

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/remuxed_read.png?raw=true)

This process allows for proper organized demultiplexing of dual-indexed outputs from single end sequencers (Ion Torrent,
PacBio, Nanopore) in QIIME2.

## Usage

### Input

The remultiplexing process can start from either a raw unmapped bam file such as what is given by an Ion Torrent
sequencer, or a fastq/fastq.gz file. Aside from the sequences, all that is needed is a .csv containing all the index tag
sequences used in the sequencing pool. This index map must follow the formatting outlined below and in the provided
example.

      sequenceFile : path to raw sequence file of type .fastq, .fastq.gz, or .bam

      indexFile    : path to .csv containing all the index tag sequences that are present in the sequencing pool. This .csv
                     should be formatted as specified below (See indexes.csv for reference)

         Column 1 - 'ID' - two digit identifiers (zero padded if single digit) for the index tags used for sequencing
         Column 2 - 'seq' - sequence of the index tag
                       *** Note that tags placed on the reverse end should be input as revers complements ***
         Column 3 - 'orientation' - F if the index was used on the Forward end, or R if on the Reverse end 

### Output

Remultiplexing will produce a single gzipped fastq containing all sequences where both a forward and reverse barcode
were found.

This format allows for immediate importing as a QIIME2 artifact of type ['MultiplexedSingleEndBarcodeInSequence']

      Return : None 

      Output : remultiplexed_seqs.fastq.gz

### Example Command Line Call

        python remultiplex.py raw_seqs.bam indexes.csv

# Index Jumping

Function: calculate the rate at which index jumps occur for a single sequencing run, and estimate the number of false
reads within the total sample pool, and for each individual sample.

Though this calculation process should be reproducible across many use cases, there are a few key requirements which
must be met for it to be accurate.

1. Dual-indexed reads (i.e. an index or barcode placed on both the forward and reverse end) which has been demultiplexed
   using QIIME2
   ![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/single_dual_index.png?raw=true)

2. Use of two or more dual-index pairs
   ![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/two_dual_index.png?raw=true)

3. 'Calibrator tags', or the use of at least one forward and one reverse index pair exclusively with each other. While
   one pair of calibraotr tags is required, we provide support for specifying multiple pairs.
   ![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/full_map_trim.png?raw=true)

## Usage

### Input

This workflow requires starting from data which has been demultiplexed utilizing QIIME2. Specifically, we start with a
data type SampleData[SequencesWithQuality], though it isn't necessary for the data to have quality values included.

    demultiplexed_seqs  : path to demultiplexed QIIME2 qza file of data type SampleData[SequencesWithQuality]

    sample_map          : path to tab delmited QIIME2 sample map file containing ALL possible index combinations
                          and a True_False indicator column (See example provided in repo)

    calibrator_tag_pairs: pairs of calibrator tags, input as either a tuple, or list of tuples, with each
                          each index being a 2 digit zero padded string.
                          ex: [('01', '11')] or [('01', '11'), ('02','12')]

### Output

    Return : Integer value of maximum number of index jumps (false reads) expected in a single sample

    Output : Expected_False_Reads_Per_Index.csv file containing number of false reads expected in EACH sample
            log.txt file containing summary statistics

### Example Command Line Call

    python index_jump.py demultiplexed_seqs.qza Sample_Map.txt 01,11

## Detailed Index Jump Calculation Workflow

The jump rate of each calibrator tag is calculated by taking the false reads with that tag, and dividing them by the
total number of that tag present in the pool.

In this instance a false read is any read which has an F01 tag and not an R11 tag, or any read with an R11 tag and not
an F01. The jump rate for the F01 tag and R11 tag are then averaged as a best estimate for the overall rate at which any
index jumps occur.

This method of calculating jump rate based off of calibrator tags offers increased accuracy when compared to calculating
jump rate based off of just the total false reads within a pool.

In using calibrator tags we are able to detect every instance in which a tag jump occurred in regard to a single index
tag (barring any jumps that don't change sample assignment). In calculating this metric using non-calibrator tags, we
lose this level of clarity. For example, if in the above example F01R12 and F01R15 were both True samples, we would not
be able to flag the jumps which occurred in these samples as false, and would conclude that no jumps occurred within our
pool. This example, while extreme, demonstrates how with each loss of a calibrator pair, our estimate of the jump rates
becomes less exact.

Ultimately this estimate is conservative, as it doesn't take into account indexes which jump but don't actually change
the sample assignment, such as the 11 index jumping from an F01R11 sample to another F01R11 sample, nor does it account
for the potential for indexes to jump at differential rates.

With this rate, we can then get an estimate for the total number of false reads in the data set by multiplying the total
read count by the jump rate. This can be used to give a good estimate of the overall percentage of true and false reads
in the data set.

Additionally, we generate a useful table (Expected_False_Read_Per_Index.csv) to assist in setting a per-sample filtering
level by listing the expected number of false reads that exist in each sample.

# PerSampleFiltering

Function: Filters reads out of a QIIME2 feature table according to a minimum read count requirement *per sample*

## Usage

### Input

      feature_table    : path to QIIME2 to feature table 

      filtering_integer: Either an integer for even filtering across samples, or path to the 
                         Expected_False_Reads_Per_Index.csv output by index_jump.py

### Output

      Return: Frequency filtered QIIME feature table of type FeatureTable[Frequency]

      Output: 'freq_filt_table.qza' QIIME artifact of type FeatureTable[Frequency]

### Example Command Line Call

    python per_sample_filtering.py feature_table.qza Expected_False_Reads_Per_Index.csv

Starting from a QIIME2 feature table consisting of 85 uniquely indexed samples. Within each sample are a number of
different features with individually recorded frequencies. Below is a histogram of a single sample with index F02R16,
which shows the frequency of each feature within the sample.

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/pre_filter.png?raw=true)

After filtering, any feature with a frequency of less than 5 are filtered out. This is carried out for each sample in
the pooled feature table.

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/post_filter.png?raw=true)
