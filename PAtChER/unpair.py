"""
Tool to remove read pairing information from SAM/BAM
"""
import sys
import argparse
import pysam


def parse_args(argv):
    """
    Parses command line options when module is run as main
    """
    args = argparse.ArgumentParser()
    args.add_argument('-i','--input', required=True,
                      metavar='<sam/bam>',
                      help='input BAM/SAM file')
    args.add_argument('-o','--output', required=True,
                      metavar='<sam/bam>',
                      help='output BAM/SAM file')
    args.add_argument('-u','--uncompressed', action="store_true",
                      help='write uncompressed (SAM) format')
    return(args.parse_args())


def main():
    """
    Removes read pairing information from a SAM/BAM file
    """
    args = parse_args(sys.argv[1:])
    input_filename = args.input
    output_filename = args.output
    sam = args.uncompressed

    input_bam = pysam.AlignmentFile(input_filename, 'r')
    if sam:
        output_bam = pysam.AlignmentFile(output_filename, 'w', header=i.header)
    else:
        output_bam = pysam.AlignmentFile(output_filename, 'wb', header=i.header)

    for segment in input_bam:
        if segment.is_reverse:
            segment.flag = 16
        else:
            segment.flag = 0
        segment.next_reference_id = -1  # set RNEXT to *
        segment.next_reference_start = -1  # set PNEXT to 0
        segment.template_length = 0  # set TLEN to 0
        output_bam.write(segment)

    input_bam.close()
    output_bam.close()
