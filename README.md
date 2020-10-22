# PAtChER

![logo](PAtChER_logo.png?raw=true "Logo")

PAtChER is a tool that uses HiChIP data to guide the assignment of multimapping reads from the ChIP step to unique locations in the genome. It thus generates a ChIP-seq output with increased genome coverage. 

## Installation

PAtChER is written in `python` and is only compatible with v3 and above. It makes use of `minimap2` to perform the mapping.

We highly recommend using virtual environments.

```bash
git clone git@github.com:MBrancoLab/PAtChER.git
virtualenv --include-lib patcherenv
. patcherenv/bin/activate
pip install ./
```

## Running

Running PAtChER is a simple then as

```bash
. patcherenv/bin/activate
patcher -h
```

Which returns

```
usage: PAtChER is a tool to help re-assign non uniquely mapping reads within a HiChIP experiment.
       [-h] -g GENOME -o OUTPUT -r1 READ1 -r2 READ2 [-d DISTANCE] [-D] [-t THREADS] [-c CUT_SITE] [-l MIN_LEN] [-b]

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Set the reference sequence (Must be FASTA)
  -o OUTPUT, --output OUTPUT
                        Set output filename
  -r1 READ1, --read1 READ1
                        Set the read sequence for 1 end (Must be FASTQ/Can be gzipped)
  -r2 READ2, --read2 READ2
                        Set the read sequence for 2 end (Must be FASTQ/Can be gzipped)
  -d DISTANCE, --distance DISTANCE
                        Set distance to map around total is +/- value
  -D, --debug           Turn on debug
  -t THREADS, --threads THREADS
                        Set the number of threads
  -c CUT_SITE, --cut-site CUT_SITE
                        Cut Site. Defaults to GATC. This will cut GATCGATC
  -l MIN_LEN, --min-len MIN_LEN
                        Minimum length read to keep for mapping. Defaults to 20
  -b, --bam             Output results in BAM format
```

## Example

A small example dataset is available in the `example_data` folder. To produce result in `SAM` format run the following command

```bash
patcher -g tests/data/test_genome.fa -r1 tests/data/ESC4_R1_sample.fq -r2 tests/data/ESC4_R2_sample.fq -o tests/data/output.sam
```

To run with multiple threads:

```bash
patcher -g tests/data/test_genome.fa -r1 tests/data/ESC4_R1_sample.fq -r2 tests/data/ESC4_R2_sample.fq -o tests/data/output.sam -t 4
```

To run output in `BAM` format:

```bash
patcher -g tests/data/test_genome.fa -r1 tests/data/ESC4_R1_sample.fq -r2 tests/data/ESC4_R2_sample.fq -o tests/data/output.bam -b
```

## Unpair alignments

From paired-end HiChIP data, PAtChER essentially produces single-end ChIP-seq alignments. Although the pairing information is preserved and could be useful, unpairing reads may be necessary for correct processing of the output in downstream applications (e.g., to generate bigwig tracks). For this purpose, the PAtChER SAM/BAM output may be processed by the unpair.py tool:

```
usage: unpair [-h] -i <sam/bam> -o <sam/bam> [-u]

optional arguments:
  -h, --help    show this help message and exit
  -i <sam/bam>  input BAM/SAM file
  -o <sam/bam>  output BAM/SAM file
  -u            write uncompressed (SAM) format
```
