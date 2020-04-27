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
pip install ./PAtChER
```

## Running

Running PAtChER is a simple then as

```bash
. patcherenv/bin/activate
patcher
```

Which returns

```
Usage: patcher -g <target.fa> -r1 <query.fq.gz> -r2 <query.fq.gz> -o <output.txt>
Options:
	-g		Set the reference sequence (Must be FASTA)
	-r1		Set the read sequence for 1 end (Must be FASTQ/Can be gzipped)
	-r2		Set the read sequence for 2 end (Must be FASTQ/Can be gzipped)
	-t		Set the number of threads
	-d		Set distance to map around total is +/- value
	-o		Set output filename
	-c		Cut Site. Default to GATC. This will cut GATCGATC
	-l		Minimum length read to keep for mapping. Defaults to 20
	-b    Output results in BAM format.
```

## Example

A small example dataset is available in the `example_data` folder. To produce result in `SAM` format run the following command

```bash
patcher -g example_data/test_genome.fa -r1 example_data/ESC4_R1_sample.fq -r2 example_data/ESC4_R2_sample.fq -o example_data/output.sam
```

To run with multiple threads:

```bash
patcher -g example_data/test_genome.fa -r1 example_data/ESC4_R1_sample.fq -r2 example_data/ESC4_R2_sample.fq -o example_data/output.sam -t 4
```

To run output in `BAM` format:

```bash
patcher -g example_data/test_genome.fa -r1 example_data/ESC4_R1_sample.fq -r2 example_data/ESC4_R2_sample.fq -o example_data/output.bam -b
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
