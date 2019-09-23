"""Main file for PAtChER. Contains main and prints help with missing required arguments.
Provides main code to run single processor version
"""
import sys
import mappy as mp
import process_reads
import alignment
import multiproc
from output import SAMBAMWriter


def getopts(argv):
    """
    Process comand line arguments
    """
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            if argv[0][1] in ["b", "D"]:
                opts[argv[0]] = True
            else:
                opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def print_help():
    """
    Method to print out command line options
    """
    print("Usage: patcher -g <target.fa> -r1 <query.fq.gz> -r2 <query.fq.gz> -o <output.txt>")
    print("Options:")
    print("\t-g\t\tSet the reference sequence (Must be FASTA)")
    print("\t-r1\t\tSet the read sequence for 1 end (Must be FASTQ/Can be gzipped)")
    print("\t-r2\t\tSet the read sequence for 2 end (Must be FASTQ/Can be gzipped)")
    print("\t-t\t\tSet the number of threads")
    print("\t-d\t\tSet distance to map around total is +/- value")
    print("\t-o\t\tSet output filename")
    print("\t-c\t\tCut Site. Defaults to GATC. This will cut GATCGATC")
    print("\t-l\t\tMinimum length read to keep for mapping. Defaults to 20")
    print("\t-b\t\tOutput results in BAM format")
    print("\t-D\t\tTurn on debug")
    exit()


def runsingle(reffile, reads1, reads2, fname, distance, cut_site, min_len, output_type):
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
            read1 = process_reads.Read(reads1.__next__())
            read1.split_read(cut_site, min_len)
            read1.qual_trim(10, 10)
            read2 = process_reads.Read(reads2.__next__())
            read2.split_read(cut_site, min_len)
            read2.qual_trim(10, 10)
            #r1 and r2 could be None at this stage
            if read1.seq and read2.seq:
                res = alignment.map_reads(reference, read1, read2, distance)
                if res:
                    sambam_output.process_output(res, read1, read2, distance)
        except StopIteration:
            break

def main(argv):
    """
    Main PAtChER function
    """
    distance = 10
    nthreads = 1
    cut_site = "GATC"
    min_len = 20
    output_type = "SAM"
    debug = False
    myargs = getopts(argv)
    if '-g' in myargs:
        reffile = myargs["-g"]
    else:
        print_help()
    if '-o' in myargs:
        fname = myargs["-o"]
    else:
        print_help()
    if '-r1' in myargs:
        r1name = myargs["-r1"]
        reads1 = mp.fastx_read(r1name)
    else:
        print_help()
    if '-r2' in myargs:
        r2name = myargs["-r2"]
        reads2 = mp.fastx_read(r2name)
    else:
        print_help()
    if '-d' in myargs:
        distance = int(myargs["-d"])
    if '-D' in myargs:
        debug = True
    if '-t' in myargs:
        nthreads = int(myargs["-t"])
    if '-c' in myargs:
        cut_site = myargs["-c"]
    if '-l' in myargs:
        min_len = int(myargs["-l"])
    if '-b' in myargs:
        output_type = "BAM"


    print(f"Using refrence:{reffile}")
    print(f"Using read 1:{r1name}")
    print(f'Using read 2:{r2name}')
    print(f'Using distance +/-:{distance}')
    print(f"Using threads:{nthreads}")
    print(f"Using cutsite:{cut_site}")
    print(f"Writing to:{fname}")
    if nthreads > 3:
        multiproc.run(reffile, reads1, reads2, fname,
                      distance, nthreads, cut_site, min_len, output_type, debug)
    else:
        if nthreads > 1:
            print("Cannot run multithreading with less than 3 threads defaulting to single")
        runsingle(reffile, reads1, reads2, fname, distance, cut_site, min_len, output_type)
    reads1.close()
    reads2.close()
    print("Run complete")
if __name__ == "__main__":
    main(sys.argv[1:])
