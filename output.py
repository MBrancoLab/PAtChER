import pysam


def get_output_var(hit):
    out = {}
    out["ctg"] = hit.ctg
    out["r_st"] = hit.r_st
    out["r_en"] = hit.r_en
    out["mapq"] = hit.mapq
    out["cigar"] = hit.cigar_str
    return out

class SAMBAMWriter():
    """docstring for ."""

    def __init__(self, filename, genome):
        self.dct_chr = {}
        header = self.create_header(genome)
        self.sambam_file = pysam.AlignmentFile(filename, "wb", header=header)

    def create_header(self, genome):
        header = { 'HD': {'VN': '1.0'},
            'SQ': [] }
        count = 0
        for k in genome.seq_names:
            header["SQ"].append({'LN': len(genome.seq(k)), 'SN': k})
            self.dct_chr[k] = count
            count += 1
        return header

    def write_read(self, read, hit1, flag, alignmentTag, hit2=None, tlen=None):
        a = pysam.AlignedSegment()
        a.query_name = read.name
        a.query_sequence=read.seq
        a.flag = flag
        a.reference_id = self.dct_chr[hit1["ctg"]]
        a.reference_start = hit1["r_st"]
        a.mapping_quality = hit1["mapq"]
        a.cigarstring = hit1["cigar"]
        if hit2:
            a.next_reference_id = self.dct_chr[hit2["ctg"]]
            a.next_reference_start = hit2["r_st"]
            a.template_length = tlen
        a.query_qualities = read.qual
        a.tags = (("PO",alignmentTag,"A"),)
        self.sambam_file.write(a)


    def processOutput(self, res, read1, read2, distance):
        #TODO handle unpaired reads
        tag=[65,129]
        hit1=get_output_var(res[0][0])
        # There is a bug in mappy that clips the mapping but doesn't provide the right cigar
        read1.seq = read1.seq[res[0][0].q_st:res[0][0].q_en]
        read1.qual = read1.qual[res[0][0].q_st:res[0][0].q_en]
        hit2=get_output_var(res[1][0])

        if res[0][1]=="u" and res[1][1]=="u":
            if hit1["ctg"] == hit2["ctg"]:
                if hit1["r_st"] < hit2["r_st"]:
                    tlen = hit2["r_en"] - hit1["r_st"]
                else:
                    tlen = hit1["r_en"] - hit2["r_st"]
                self.write_read(read1, hit1, tag[0], res[0][1], hit2, tlen)
                self.write_read(read2, hit2, tag[1], res[1][1], hit1, tlen)
        elif res[0][1] != "u":
            if hit1["r_en"] > distance:
                tlen = hit1["r_en"] - distance
            else:
                tlen = distance - hit1["r_st"]
            hit1["r_st"] = hit2["r_st"] - distance + hit1["r_st"]
            hit1["ctg"] = hit2["ctg"]
            self.write_read(read1, hit1, tag[0], res[0][1], hit2, tlen)
            self.write_read(read2, hit2, tag[1], res[1][1], hit1, tlen)
        elif res[1][1]!="u":
            if hit2["r_en"] > distance:
                tlen = hit2["r_en"] - distance
            else:
                tlen = distance - hit2["r_st"]
            hit2["r_st"] = hit1["r_st"] - distance + hit2["r_st"]
            hit2["ctg"] = hit1["ctg"]
            self.write_read(read1, hit1, tag[0], res[0][1], hit2, tlen)
            self.write_read(read2, hit2, tag[1], res[1][1], hit1, tlen)

    def __del__(self):
        self.sambam_file.close()
