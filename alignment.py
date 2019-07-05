import mappy as mp
import math
import random

#Function to map a read par consisting of a 2D List
# a -> mappy reference
# distance -> +/- distance to search for non unique hit
def map(read1,read2,a,distance):
    hits=[[],[]]
    hits[0]=list(a.map(read1.seq))
    hits[1]=list(a.map(read2.seq))
    if len(hits[0])==1 and len(hits[1])==1:
        hits[0].append("u")
        hits[1].append("u")
        return hits
    else:
        return(process_unique_one(hits,a,read1,read2,distance))

#Function to process a unique hit
def process_unique_one(hits,a,read1,read2,distance):
    if len(hits[0])==1 and len(hits[1])>0:
        indx = 0
        seq = read2.seq
    elif len(hits[1])==1 and len(hits[0])>0:
        indx = 1
        seq = read1.seq
    else:
        return None
    #Get reference sequence for unique hit
    refseq=a.seq(hits[indx][0].ctg,hits[indx][0].r_st-distance,hits[indx][0].r_en+distance)
    if refseq:
        b= mp.Aligner(seq=refseq,preset="sr",n_threads=1)
        newhits=[]
        for hit in b.map(seq): # traverse alignments
            if hit.mlen/(len(seq)*1.0)>0.8:
                newhits.append(hit)
        #Need to fix hit.ctg and hit.r_st???
        if len(newhits)==1:
            if indx==0:
                hits[0].append("u")
                hits[1]=[newhits[0],"r"]
            else:
                hits[1].append("u")
                hits[0]=[newhits[0],"r"]
            return(hits)
        elif len(newhits)>0:
            ds=[]
            for k in newhits:
                if k.r_st>=distance:
                    ds.append(k.r_st-distance)
                else:
                    ds.append(distance-k.r_en)
            s1=0
            s2=[]
            for k in ds:
                a=math.exp(-0.8*k/50-0.6618)
                if a<0:
                    a=0
                s1+=a
                s2.append(s1)
                p=random.random()*s1
            for k in range(len(s2)):
                if p<=s2[k]:
                    hits[indx]=[newhits[k],"p"]
                    return(hits)
    return(None)
