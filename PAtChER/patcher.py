"""Main file for PAtChER. Contains main and prints help with missing required arguments.
Provides main code to run single processor version
"""
import sys
import os.path
from argparse import ArgumentParser
import mappy as mp
from .process_reads import Read
from .alignment import map_reads
from .multiproc import run
from .output import SAMBAMWriter


def runsingle(reffile, reads1, reads2, fname, distance, cut_site, min_rec, min_len, output_type):
    """Run single processor version of PAtChER"""
    print("Loading in Reference")
    reference = mp.Aligner(reffile, preset="sr")
    if not reference:
        raise Exception("ERROR: failed to load/build index file")
    print("Done.")
    sambam_output = SAMBAMWriter(fname, reference, output_type)
    print("Running Alignment")
    while True:
        try:
            read1 = Read(reads1.__next__())
            read1.split_read(cut_site, min_rec, min_len)
            read1.qual_trim(10, 10)
            read2 = Read(reads2.__next__())
            read2.split_read(cut_site, min_rec, min_len)
            read2.qual_trim(10, 10)
            #r1 and r2 could be None at this stage
            if read1.seq and read2.seq:
                res = map_reads(reference, read1, read2, distance)
                if res:
                    sambam_output.process_output(res, read1, read2, distance)
        except StopIteration:
            break

def main():
    """
    Main PAtChER function
    """
    parser = ArgumentParser(
        "PAtChER is a tool to help re-assign non uniquely mapping reads within a HiChIP experiment.")
    parser.add_argument('-g', '--genome', required=True, help="Set the reference sequence (Must be FASTA)")
    parser.add_argument('-o', '--output', required=True, help="Set output filename")
    parser.add_argument('-r1',
                        '--read1',
                        required=True,
                        help="Set the read sequence for 1 end (Must be FASTQ; can be gzipped)")
    parser.add_argument('-r2',
                        '--read2',
                        required=True,
                        help="Set the read sequence for 2 end (Must be FASTQ; can be gzipped)")
    parser.add_argument('-d',
                        '--distance',
                        type=int,
                        required=False,
                        default=10000,
                        help="Set distance from unique ends to map non-unique end. Defaults to 10kb.")
    parser.add_argument('-D', '--debug', required=False, action='store_true', help="Turn on debug")
    parser.add_argument('-t', '--threads', type=int, required=False, default=1, help="Set the number of threads")
    parser.add_argument('-c',
                        '--cut-site',
                        required=False,
                        default="GATCGATC",
                        help="Cut Site. Defaults to GATCGATC, i.e., a DpnII site after digestion, fill-in and ligation")
    parser.add_argument('-m',
                        '--min-rec',
                        type=int,
                        required=False,
                        default=5,
                        help="Minimum length of cut site to recognise when at the 3' end. Defaults to 5")
    parser.add_argument('-l',
                        '--min-len',
                        type=int,
                        required=False,
                        default=40,
                        help="Minimum read length to keep for mapping. Defaults to 40")
    parser.add_argument('-b', '--bam', required=False, action='store_true', help="Output results in BAM format")
    args = parser.parse_args()

    if os.path.exists(args.read1):
        reads1 = mp.fastx_read(args.read1)
    else:
        raise Exception(f"ERROR: could not open file {args.read1}")
    if os.path.exists(args.read2):
        reads2 = mp.fastx_read(args.read2)
    else:
        raise Exception(f"ERROR: could not open file {args.read2}")

    print(f"Using refrence: {args.genome}")
    print(f"Using read 1: {args.read1}")
    print(f'Using read 2: {args.read2}')
    print(f'Using distance +/-: {args.distance}')
    print(f"Using threads: {args.threads}")
    print(f"Using cutsite: {args.cut_site}")
    print(f"Writing to: {args.output}")
    if args.threads >= 3:
        run(args.genome, reads1, reads2, args.output,
                      args.distance, args.threads, args.cut_site, args.min_rec, args.min_len, args.bam, args.debug)
    else:
        if args.threads > 1:
            print("Cannot run multithreading with less than 3 threads. Defaulting to single")
        runsingle(args.genome, reads1, reads2, args.output, args.distance, args.cut_site, args.min_rec, args.min_len, args.bam)
    reads1.close()
    reads2.close()
    print("Run complete")
