# PAtChER

![logo](PAtChER_logo.png?raw=true "Logo")

PAtChER is a tool to help re-assign non uniquely mapping reads within a HiChIP experiment.

## Installation

PAtChER is written in `python` and is only compatible with v3 and above. It makes use of `minimap2` to perform the mapping.

We highly recommend using virtual environments.

```bash
virtualenv --include-lib patcherenv
. patcherenv/bin/activate
pip install mappy
pip install pysam
```

## Running

Running PAtChER is a simple then as

```bash
. patcherenv/bin/activate
python3 patcher.py
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
python3 patcher.py -g example_data/test_genome.fa -r1 example_data/ESC4_R1_sample.fq -r2 example_data/ESC4_R2_sample.fq -o example_data/output.sam
```

To run with multiple threads:

```bash
python3 patcher.py -g example_data/test_genome.fa -r1 example_data/ESC4_R1_sample.fq -r2 example_data/ESC4_R2_sample.fq -o example_data/output.sam -t 4
```

To run output in `BAM` format:

```bash
python3 patcher.py -g example_data/test_genome.fa -r1 example_data/ESC4_R1_sample.fq -r2 example_data/ESC4_R2_sample.fq -o example_data/output.bam -b
```
