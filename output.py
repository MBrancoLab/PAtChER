"""
Module to handle output to SAM or BAM
"""
import pysam


class SAMBAMWriter():
    """Class for writing output to SAM or BAM format"""

    def __init__(self, filename, genome, output_type):
        self.dct_chr = {}
        header = self.create_header(genome)
        if output_type == "SAM":
            self.sambam_file = pysam.AlignmentFile(filename, "w", header=header)
        elif output_type == "BAM":
            self.sambam_file = pysam.AlignmentFile(filename, "wb", header=header)
        else:
            raise Exception(f"Unknown output_type: {output_type}. Should be SAM or BAM")

    def create_header(self, genome):
        """
        Method to produce the header for SAM or BAM file
        """
        header = {'HD': {'VN': '1.0'},
                  'SQ': []}
        count = 0
        for k in genome.seq_names:
            header["SQ"].append({'LN': len(genome.seq(k)), 'SN': k})
            self.dct_chr[k] = count
            count += 1
        return header

    def write_read(self, read, hit1, flag, alignment_tag1, alignment_tag2, alignment_tag3, alignment_tag4,
                   tlen=None, hit2=None):
        """
        Method which writes data to file
        """
        aligned_segment = pysam.AlignedSegment()
        aligned_segment.query_name = read.name
        aligned_segment.query_sequence = read.seq
        aligned_segment.flag = flag
        aligned_segment.reference_id = self.dct_chr[hit1["ctg"]]
        aligned_segment.reference_start = hit1["r_st"]
        aligned_segment.mapping_quality = hit1["mapq"]
        aligned_segment.cigarstring = hit1["cigar"]

        if hit2:
            aligned_segment.next_reference_id = self.dct_chr[hit2["ctg"]]
            aligned_segment.next_reference_start = hit2["r_st"]
            aligned_segment.template_length = tlen
        aligned_segment.query_qualities = read.qual
        aligned_segment.tags = (("XM", alignment_tag1, "Z"), ("XR", alignment_tag2, "Z"), ("XG", alignment_tag3, "Z"),
                                ("PO", alignment_tag4, "A"),)
        self.sambam_file.write(aligned_segment)

    def process_paired_alignment(self, res, read1, read2):
        """
        Method to process a paired alignment
        """
        tlen = 0
        flag = [65, 129]
        if res[0][0]["strand"] == -1:
            flag[0] += 16
            flag[1] += 32
        if res[1][0]["strand"] == -1:
            flag[0] += 32
            flag[1] += 16

        if res[0][0]["strand"] == 1:
            tag3_1 = "CT"
        if res[0][0]["strand"] == -1:
            tag3_1 = "GA"
        if res[1][0]["strand"] == 1:
            tag3_2 = "GA"
        if res[1][0]["strand"] == -1:
            tag3_2 = "CT"

        # There is a bug in mappy that clips the mapping but doesn't provide the right cigar
        read1.seq = read1.seq[res[0][0]["q_st"]:res[0][0]["q_en"]]
        read1.qual = read1.qual[res[0][0]["q_st"]:res[0][0]["q_en"]]
        read2.seq = read2.seq[res[1][0]["q_st"]:res[1][0]["q_en"]]
        read2.qual = read2.qual[res[1][0]["q_st"]:res[1][0]["q_en"]]
        hit1 = res[0][0]
        hit2 = res[1][0]
        if hit1["ctg"] == hit2["ctg"]:
            if hit1["r_st"] < hit2["r_st"]:
                tlen = hit2["r_en"] - hit1["r_st"]
            else:
                tlen = hit1["r_en"] - hit2["r_st"]

        self.write_read(read1, hit1, flag[0], res[0][2], "CT", tag3_1, res[0][1], tlen, hit2)
        self.write_read(read2, hit2, flag[1], res[1][2], "GA", tag3_2, res[1][1], tlen, hit1)

    def process_output(self, res, read1, read2):
        """
        Take output from alignment (res) and reads (read1/read2)
        and process them to SAM/BAM
        """

        flag = 0
        if len(res[1]) == 0:
            hit = res[0][0]
            read1.seq = read1.seq[res[0][0]["q_st"]:res[0][0]["q_en"]]
            read1.qual = read1.qual[res[0][0]["q_st"]:res[0][0]["q_en"]]
            tag2 = "CT"
            tag3 = "CT"
            if res[0][0]["strand"] == -1:
                flag += 16
                tag3 = "GA"
            self.write_read(read1, hit, flag, res[0][2], tag2, tag3, res[0][1])

        elif len(res[0]) == 0:
            hit = res[1][0]
            read2.seq = read2.seq[res[1][0]["q_st"]:res[1][0]["q_en"]]
            read2.qual = read2.qual[res[1][0]["q_st"]:res[1][0]["q_en"]]
            tag2 = "GA"
            tag3 = "GA"
            if res[1][0]["strand"] == -1:
                flag += 16
                tag3 = "CT"
            self.write_read(read2, hit, flag, res[1][2], tag2, tag3, res[1][1])

        elif len(res[0]) == 3 and len(res[1]) == 3:
            self.process_paired_alignment(res, read1, read2)

        else:
            raise Exception(f"Incorrect number of alignment hits: {len(res[0])} and {len(res[1])}")

    def close_file(self):
        self.sambam_file.close()