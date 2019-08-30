"""
Module containing the class for Read
"""
import re

class read:
    def __init__(self, dt):
        self.name = dt[0]
        self.seq = dt[1]
        self.qual_str = dt[2]
        self.qual = [ord(x) - 33 for x in dt[2]]
        self.CT = self.seq_conversion("C", "T") # C-to-T version read
        self.GA = self.seq_conversion("G", "A") # G-to-A version read

    def seq_conversion(self, base1, base2): # bisulfite convertion of reads
        seq = self.seq.replace(base1, base2)
        return (seq.replace(base1.lower(), base2.lower()))

    def splitRead(self, recseq, minlen):
        hitsite = -1
        for hit in re.finditer(f"{recseq}{recseq[0]}", self.seq):
            hit_pos = hit.start()
            if hit_pos != -1:
                site = True
                for counter, value in enumerate(recseq):
                    if hit_pos + counter + len(recseq) < len(self.seq):
                        if self.seq[hit_pos + counter + len(recseq)] != value:
                            site = False
                if site:
                    hitsite = hit_pos
                    break
        if hitsite != -1:
            if hitsite < minlen:
                self.seq = None
                self.qual_str = None
                self.qual = None
                self.CT = None
                self.GA = None
            else:
                self.seq = self.seq[:(hitsite + len(recseq))]
                self.qual_str = self.qual_str[:(hitsite + len(recseq))]
                self.qual = self.qual[:(hitsite + len(recseq))]
                self.CT = self.CT[:(hitsite + len(recseq))]
                self.GA = self.GA[:(hitsite + len(recseq))]