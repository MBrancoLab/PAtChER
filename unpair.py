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
    args.add_argument('-i', required=True,
                      metavar='<sam/bam>',
                      help='input BAM/SAM file')
    args.add_argument('-o', required=True,
                      metavar='<sam/bam>',
                      help='output BAM/SAM file')
    args.add_argument('-u', action="store_true",
                      help='write uncompressed (SAM) format')
    return(args.parse_args())
        

def unpair(ifile, ofile, sam):
    """
    Removes read pairing information from a SAM/BAM file
    """
    i = pysam.AlignmentFile(ifile, 'r')
    if sam:
        o = pysam.AlignmentFile(ofile, 'w', header = i.header)
    else:
        o = pysam.AlignmentFile(ofile, 'wb', header = i.header)

    for segment in i:
        if segment.is_reverse:
            segment.flag = 16
        else:
            segment.flag = 0
        segment.next_reference_id = -1  # set RNEXT to *
        segment.next_reference_start = -1  # set PNEXT to 0
        segment.template_length = 0  # set TLEN to 0
        o.write(segment)

    i.close()
    o.close()


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    unpair(args.i, args.o, args.u)
