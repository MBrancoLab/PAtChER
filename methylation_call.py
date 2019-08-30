import mappy as mp
import re

def call(res, read, reference, readN):
    """
        Method to process methylation calling
    """
    read_seq = read.seq[res[0]["q_st"]:res[0]["q_en"]]
    ref_seq = reference.seq(res[0]['ctg'], res[0]['r_st'], res[0]['r_en'])
    call = '.' * len(read_seq)
    call = list(call)
    sites = []
    if readN == 1:  # read1
        if res[0]['strand'] == 1: # when read1 is on the forward strand
            if len(read_seq) < len(ref_seq): # keep same length of read_seq and ref_seq
                ref_seq = ref_seq[:len(read_seq)]
            else:
                read_seq = read_seq[:len(ref_seq)]

            for hit in re.finditer("CG", ref_seq):
                site = hit.start()
                if site != -1:
                    sites.append(site)

            for x in sites:
                if read_seq[x] == 'C':
                    call[x] = 'Z' # methylated C
                elif read_seq[x] == 'T':
                    call[x] = 'z' # unmethylated C
                else:
                    call[x] = '.'

        elif res[0]['strand'] == -1: # when read1 is on the reverse strand
            ref_seq = mp.revcomp(ref_seq) # comparing read_seq with complementary ref_seq
            if len(read_seq) < len(ref_seq):
                ref_seq = ref_seq[:len(read_seq)]
            else:
                read_seq = read_seq[:len(ref_seq)]

            for hit in re.finditer("GC", ref_seq):
                site = hit.start()
                if site != -1:
                    sites.append(site)

            for x in sites:
                if read_seq[x + 1] == 'C':
                    call[x + 1] = 'Z'
                elif read_seq[x + 1] == 'T':
                    call[x + 1] = 'z'
                else:
                    call[x + 1] = '.'

    else:  # read2
        if res[0]['strand'] == 1: # when read2 is on the forward strand
            if len(read_seq) < len(ref_seq):
                ref_seq = ref_seq[:len(read_seq)]
            else:
                read_seq = read_seq[:len(ref_seq)]

            for hit in re.finditer("GC", ref_seq):
                site = hit.start()
                if site != -1:
                    sites.append(site)

            for x in sites:
                if read_seq[x] == 'G':
                    call[x] = 'Z'
                elif read_seq[x] == 'A':
                    call[x] = 'z'
                else:
                    call[x] = '.'

        elif res[0]['strand'] == -1: # when read2 is on the reverse strand
            ref_seq = mp.revcomp(ref_seq)

            if len(read_seq) < len(ref_seq):
                ref_seq = ref_seq[:len(read_seq)]
            else:
                read_seq = read_seq[:len(ref_seq)]

            for hit in re.finditer("CG", ref_seq):
                site = hit.start()
                if site != -1:
                    sites.append(site)

            for x in sites:
                if read_seq[x + 1] == 'G':
                    call[x + 1] = 'Z'
                elif read_seq[x + 1] == 'A':
                    call[x + 1] = 'z'
                else:
                    call[x + 1] = '.'

    call = "".join(call)
    return (call)