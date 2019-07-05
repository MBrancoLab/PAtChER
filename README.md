# PAtChER

![logo](patcher_logo.png?raw=true "Logo")

PAtChER is a tool to help re-assign non uniquely mapping reads within a HiChIP experiment.

## Installation

PAtChER is written in `python` and is only compatiable with v3 and above. It makes use of `minimap2` to perform the mapping.

We highly recommend using virtual enviroments.

```bash
virtualenv --include-lib patcherenv
. patcherenv/bin/activate
pip install mappy
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
	-c		Cut Site. Defaults to GATC. This will cut GATCGATC
	-l		Minumum length read to keep for mapping. Defaults to 20
```

**TODO ADD EXAMPLE WALKTHROUGH WITH DATA**