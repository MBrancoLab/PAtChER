"""
Module containing the functions to align reads to the genome and process non unique alignments
"""
import math
import random
import mappy as mp

def get_output_var(hit):
    """
    Method to convert hit to dictonary to make mutable
    """
    out = {}
    out["ctg"] = hit.ctg
    out["r_st"] = hit.r_st
    out["r_en"] = hit.r_en
    out["q_st"] = hit.q_st
    out["q_en"] = hit.q_en
    out["mapq"] = hit.mapq
    out["cigar"] = hit.cigar_str
    out["strand"] = hit.strand
    return out


#Function to map a read par consisting of a 2D List
# a -> mappy reference
# distance -> +/- distance to search for non unique hit
def map_reads(reference, read1, read2, distance):
    """
    Method to map read pair
    reference -> mappy reference
    distance -> to search around unique mapped read with non unique pair
    """
    hits = [[], []]
    hits[0] = list(reference.map(read1.seq))
    hits[1] = list(reference.map(read2.seq))
    if len(hits[0]) == 1 and len(hits[1]) == 1:
        out = [[get_output_var(hits[0][0]), "u"], [get_output_var(hits[1][0]), "u"]]
    elif len(hits[0]) == 1 and len(hits[1]) == 0:
        out = [[get_output_var(hits[0][0])], []]
    elif len(hits[0]) == 0 and len(hits[1]) == 1:
        out = [[], [get_output_var(hits[1][0])]]
    elif len(hits[0]) == 1 and len(hits[1]) > 0:
        out = process_unique_one(reference, hits, read1, read2, distance)
    elif len(hits[1]) == 1 and len(hits[0]) > 0:
        out = process_unique_one(reference, hits, read1, read2, distance)
    else:
        out = None

    return out

def process_unique_one(reference, hits, read1, read2, distance):
    """
    Method to process non unique mapping hits
    """
    if len(hits[0]) == 1 and len(hits[1]) > 0:
        indx = 0
        seq = read2.seq
        out = [[get_output_var(hits[0][0]), "u"], []]
    elif len(hits[1]) == 1 and len(hits[0]) > 0:
        indx = 1
        seq = read1.seq
        out = [[], [get_output_var(hits[1][0]), "u"]]
    #Get reference sequence for unique hit
    refseq = reference.seq(hits[indx][0].ctg,
                           hits[indx][0].r_st - distance, hits[indx][0].r_en + distance)
    if refseq:
        local_reference = mp.Aligner(seq=refseq, preset="sr", n_threads=1)
        new_hits = []
        for hit in local_reference.map(seq): # traverse alignments
            if hit.mlen/(len(seq) * 1.0) > 0.8:
                new_hits.append(hit)
        # Need to fix hit.ctg and hit.r_st???
        if len(new_hits) == 1:
            if indx == 0:
                out[1] = [get_output_var(new_hits[0]), "r"]
                out[1][0]["ctg"] = out[0][0]["ctg"]
            else:
                out[0] = [get_output_var(new_hits[0]), "r"]
                out[0][0]["ctg"] = out[1][0]["ctg"]
        elif len(new_hits) > 0:
            distance_list = []
            for new_hit in new_hits:
                if new_hit.r_st >= distance:
                    distance_list.append(new_hit.r_st - distance)
                else:
                    distance_list.append(distance - new_hit.r_en)
            probability_list = []
            for dist in distance_list:
                scale_probability = math.exp(-0.8 * dist/50 - 0.6618)
                probability_list.append(scale_probability)
            pindex = random.choices(range(len(probability_list)), weights=probability_list)[0]
            if indx == 0:
                out[1] = [get_output_var(new_hits[pindex]), "p"]
                out[1][0]["ctg"] = out[0][0]["ctg"]
            else:
                out[0] = [get_output_var(new_hits[pindex]), "p"]
                out[0][0]["ctg"] = out[1][0]["ctg"]
    return out
