# MetaPlex 

*BROADLY DESCRIBE THE METAPLEX WORKFLOW*

# Remultiplexing 

Including unique indexes (or barcodes) on individual samples in a pooled sequencing run is a great way to optimize
efficiency when performing next-generation sequencing. Proper downstream analysis of these sequencing runs is dependent
on properly assigning each sequence to its source sample through organized demultiplexing.

Here we provide a tool which allows for proper organized demultiplexing of dual-indexed outputs from single end
sequencers (Ion Torrent, PacBio, Nanopore) in QIIME2 through a 'remultiplexing' process.

In this example, 'Ion Xpress Barcodes' were utilized to uniquely index each sample on both the forward and the reverse
end, with indexes 1 through 10 being placed on the Forward end, and 11 through 20 being placed on the Reverse end. In
order to easily process these reads in QIIME2, both of these indexes must be linked, and placed together at the front of
the reads.

## Dependencies

Command line tools:

- samtools (if starting from .bam file)
- cutadapt

Python3 with the following packages installed:

- subprocess
- os
- sys
- pandas

## Usage

The remultiplexing process can start from either a raw unmapped bam file such as what is given by an Ion Torrent
sequencer, or a fastq/fastq.gz file. Aside from the sequences, all that is needed is a .csv containing all the index tag
sequences used in the sequencing pool. This index map must follow the formatting of the one provided in this repo.

    Column 1: ID - two digit identifiers (zero padded if single digit) for the index tags used for sequencing
    Column 2: seq - sequence of the index tag
                    *** Note that tags placed on the reverse end should be input as revers complements ***
    Column 3: orientation - F if the index was used on the Forward end, or R if on the Reverse end 

### In a Jupyter Notebook

To import the function for usage in a Jupyter notebook:

        from remultiplexing import remultiplex

### As a command line tool

To call from the command line:

        python3 remultiplex.py *SEQUENCE_FILE* *INDEX_CSV*

### Output

Running this will produce a single gzipped fastq containing all sequences where a forward and reverse barcode was found.

This file be output as ***remultiplexed_seqs.fastq.gz***

This format allows for immediate importing as a QIIME2 artifact of type ['MultiplexedSingleEndBarcodeInSequence']

# Index Jumping 

This tool is used to calculate the rate at which index jumps occur for a single sequencing run, and uses this to 
estimate the number of false reads within a sample pool. Though this calculation process should be reproducible across 
many use cases, there are a few key requirements which must be met for it to be accurate.

1. Dual-indexed reads (i.e. an index or barcode placed on both the forward and reverse end) which has been 
demultiplexed using QIIME2
![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/single_dual_index.png?raw=true)

2. Use of two or more dual-index pairs 
![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/two_dual_index.png?raw=true)

3. 'Calibrator tags', or the use of at least one forward and one reverse index pair exclusively with each other.
While one pair of calibraotr tags is required, we provide support for specifying  multiple pairs.
![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/full_map_trim.png?raw=true)

Provided in this repo is a python script which can be called from the command line or used in a Jupyter Notebook.
A stand alone Jupyter Notebook outlining the basic process is also available. 

## Dependencies

Script should be used in an active QIIME2 (>= 2021.8) conda environment 

Additional python imports:
* pathlib 
* sys
* shutil
* tempfile
* os
* math
* pandas


## Inputs

This workflow requires starting from data which has been demultiplexed utilizing QIIME2. Specifically, we start with a 
data type SampleData[SequencesWithQuality], though it isn't necessary for the data to have quality values included.

    demultiplexed_seqs  : path to demultiplexed QIIME2 qza file of data type SampleData[SequencesWithQuality]

    sample_map          : path to tab delmited QIIME2 sample map file containing ALL possible index combinations
                          and a True_False indicator column (See example provided in repo)

    calibrator_tag_pairs: pairs of calibrator tags, input as either a tuple, or list of tuples, with each
                          each index being a 2 digit zero padded string.
                          ex: [('01', '11')] or [('01', '11'), ('02','12')]


## Outputs

    Return: Integer value of maximum number of index jumps (false reads) expected in a single sample

    Output: Expected_False_Reads_Per_Index.csv file containing number of false reads expected in EACH sample
            log.txt file containing summary statistics

## Example Command Line Call 

    python index_jump.py demultiplexed_seqs.qza Sample_Map.txt 01,11


## Detailed Workflow

The jump rate of each calibrator tag is calculated by taking the false reads with that tag, and dividing them by the 
total number of that tag present in the pool. 

In this instance a false read is any read which has an F01 tag and not an R11 tag, or any read with an R11 tag and not 
an F01. The jump rate for the F01 tag and R11 tag are then averaged as a best estimate for the overall rate at which 
any index jumps occur.

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

Additionally, we generate a useful table (Expected_False_Read_Per_Index.csv) to assist in setting a per-sample 
filtering level by listing the expected number of false reads that exist in each sample. 

For more detailed examples using a real sequencing output, consult the Jupyter-Notebook provided in this repo.


# PerSampleFiltering

Here is a short script and jupyter-notebook which walks through the process of filtering reads/features out of a QIIME2 feature table according to a minimum read count requirement. 

In some instances it may be necessary to filter out any read per-sample which doesn't appear at a minimum threshold count. Currently QIIME2 offers no native way of performing this action. 

In this example I used a small subset of COI barcode sequences from a previous sequencing project to demonstrate filtering out features in EACH sample which don't appear more than 20 times. 

Here, we start with a qiime2 feature table consisting of ~70 uniquely indexed samples. Within each sample are a number of different features with indivudally recorded frequencies. Below is a histogram of a single sample with index F02R16, which shows the frequency of occurence of each feature within the sample. 

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/pre_filter.png?raw=true)

After we perform our filtering steps outlined in the jupyter-notebook, we can see below that any feature with a frequency of occurence less than or equal to 20 are filtered out. This is carried out for each sample in the pooled feature table. 

![alt text](https://github.com/NGabry/MetaPlex/blob/main/images/post_filter.png?raw=true)

Provided in this repo are the feature-table.qza, as well as a feature-table.biom file, to demonstrate a couple of different starting points for this potential workflow

## Using the .py script 

Additinoally I have included a script which can be executed from the command line for faster processing.    

To run this script the following dependencies must be installed within an active Qiime2 environment
* biom
* pandas

When runnign the script, just specify the integer value at which you'd like to filter and state this as a command line argument! 
For example, to filter out any feature which appears less than 20 times, type the following:

    python per_sample_filtering.py feature_table.qza 20
    
