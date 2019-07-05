import re

class read:
    #Mappy produces a tuple for each read
    def __init__(self,dt):
        self.name = dt[0]
        self.seq = dt[1]
        self.qual_str = dt[2]
        self.qual = [ord(x)-33 for x in dt[2]]

    def splitRead(self,recseq,minlen):
        hitsite= -1
        for hit in re.finditer(f"{recseq}{recseq[0]}", self.seq):
            p = hit.start()
            if p!=-1:
                site=True
                for k in range(len(recseq)):
                    if p+k+len(recseq) < len(self.seq):
                        if self.seq[p+k+len(recseq)]!=recseq[k]:
                            site=False
                if site:
                    hitsite=p
                    break
        if hitsite!=-1:
            if hitsite<minlen:
                self.seq=None
                self.qual=None
            else:
                self.seq=self.seq[:(hitsite+len(recseq))]
                self.qual=self.qual[:(hitsite+len(recseq))]

    def qualtrim(self,cutoff_front,cutoff_back):
        #Taken from cutadapt which in tern was taken from bwa_trim_read
        #Check that we havent already lost the read e.g. qual is None
        if self.qual:
            start=0
            stop=len(self.qual)
            s=0
            max_qual=0
            for i in range(len(self.qual)):
                s+= cutoff_front-self.qual[i]
                if s<0:
                    break
                if s>max_qual:
                    max_qual = s
                    start = i+1

            s=0
            max_qual=0
            for i in reversed(range(len(self.qual))):
                s+=cutoff_front-self.qual[i]
                if s<0:
                    break
                if s>max_qual:
                    max_qual = s
                    stop = i
            if start>=stop:
                self.qual = None
                self.seq = None
            else:
                self.seq=self.seq[start:stop]
                self.qual=self.qual[start:stop]
                self.qual_str=self.qual_str[start:stop]


                
                
                