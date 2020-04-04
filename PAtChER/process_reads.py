"""
Module containing the class for Read
"""
import re

class Read:
    """
    Class to store reads read in using the mappy function.
    Provide methods for quality trimming and spliting of reads
    based on a recognition sequence
    """
    def __init__(self, dt):
        self.name = dt[0]
        self.seq = dt[1]
        self.qual_str = dt[2]
        self.qual = [ord(x)-33 for x in dt[2]]

    def split_read(self, recseq, minlen):
        """
        Method to split reads based on a recognition sequence
        """
        hitsite = -1
        for hit in re.finditer(f"{recseq}{recseq[0]}", self.seq):
            hit_pos = hit.start()
            if hit_pos != -1:
                site = True
                for counter, value in enumerate(recseq):
                    if hit_pos+counter+len(recseq) < len(self.seq):
                        if self.seq[hit_pos+counter+len(recseq)] != value:
                            site = False
                if site:
                    hitsite = hit_pos
                    break
        if hitsite != -1:
            if hitsite < minlen:
                self.seq = None
                self.qual = None
            else:
                self.seq = self.seq[:(hitsite+len(recseq))]
                self.qual = self.qual[:(hitsite+len(recseq))]

    def qual_trim(self, cutoff_front, cutoff_back):
        """
        Method to perform quality trimming on reads.
        Taken from cutadapt which in tern was taken from bwa_trim_read
        """
        # Check that we havent already lost the read e.g. qual is None
        if self.qual:
            start = 0
            stop = len(self.qual)
            current_qual = 0
            max_qual = 0
            for i in range(len(self.qual)):
                current_qual += cutoff_front-self.qual[i]
                if current_qual < 0:
                    break
                if current_qual > max_qual:
                    max_qual = current_qual
                    start = i+1

            current_qual = 0
            max_qual = 0
            for i in reversed(range(len(self.qual))):
                current_qual += cutoff_back-self.qual[i]
                if current_qual < 0:
                    break
                if current_qual > max_qual:
                    max_qual = current_qual
                    stop = i
            if start >= stop:
                self.qual = None
                self.seq = None
            else:
                self.seq = self.seq[start:stop]
                self.qual = self.qual[start:stop]
                self.qual_str = self.qual_str[start:stop]
