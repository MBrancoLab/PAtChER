"""
Module containing the Alignment for paired-end read
"""

import mappy as mp
import math
import random


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


def map_paired(read1, read2, a, top_a, bottom_a, distance):
    hits = [[], [], [], []]

    for hit1 in top_a.map(read1.CT):  # r1.CT against top_ref
        if hit1.mlen / (len(read1.seq) * 1.0) > 0.8:
            hits[0].append(hit1)

    for hit2 in bottom_a.map(read1.CT):  # r1.CT against bottom_ref
        if hit2.mlen / (len(read1.seq) * 1.0) > 0.8:
            hits[1].append(hit2)

    for hit3 in top_a.map(read2.GA):  # r2.GA against top_ref
        if hit3.mlen / (len(read2.seq) * 1.0) > 0.8:
            hits[2].append(hit3)

    for hit4 in bottom_a.map(read2.GA):  # r2.GA against bottom_ref
        if hit4.mlen / (len(read2.seq) * 1.0) > 0.8:
            hits[3].append(hit4)

    if len(hits[0]) == len(hits[2]) == 1: # Unique hit: r1.CT against top_ref & r2.GA against top_ref
        out = [[get_output_var(hits[0][0]), 'u'], [get_output_var(hits[2][0]), 'u']]

    elif len(hits[0]) == len(hits[3]) == 1: # Unique hit: r1.CT against top_ref & r2.GA against bottom_ref
        out = [[get_output_var(hits[0][0]), 'u'], [get_output_var(hits[3][0]), 'u']]

    elif len(hits[1]) == len(hits[3]) == 1: # Unique hit: r1.CT against bottom_ref & r2.GA against bottom_ref
        out = [[get_output_var(hits[1][0]), 'u'], [get_output_var(hits[3][0]), 'u']]

    elif len(hits[1]) == len(hits[2]) == 1: # Unique hit: r1.CT against bottom_ref & r2.GA against top_ref
        out = [[get_output_var(hits[1][0]), 'u'], [get_output_var(hits[2][0]), 'u']]

    elif len(hits[0]) == 1 and len(hits[1]) == len(hits[2]) == len(hits[3]) == 0: # single-end hit: r1.CT against top_ref
        out = [[get_output_var(hits[0][0]), 's'], []]

    elif len(hits[1]) == 1 and len(hits[0]) == len(hits[2]) == len(hits[3]) == 0: # single-end hit: r1.CT against bottom_ref
        out = [[get_output_var(hits[1][0]), 's'], []]

    elif len(hits[2]) == 1 and len(hits[0]) == len(hits[1]) == len(hits[3]) == 0: # single-end hit: r2.GA against top_ref
        out = [[], [get_output_var(hits[2][0]), 's']]

    elif len(hits[3]) == 1 and len(hits[0]) == len(hits[1]) == len(hits[2]) == 0: # single-end hit: r2.GA against bottom_ref
        out = [[], [get_output_var(hits[3][0]), 's']]

    elif len(hits[0]) == 1 and len(hits[2]) > 1 and len(hits[3]) > 1: # r1.CT against top_ref unique hit, r2 non-unique hit
        h1 = [hits[0], hits[2]]
        h2 = [hits[0], hits[3]]
        out1 = process_unique_one(h1, read1, read2, a, top_a, distance) # unique r1 hit and multiple r2 against top_ref to process_unique_one
        out2 = process_unique_one(h2, read1, read2, a, bottom_a, distance) # unique r1 hit and multiple r2 against bottom_ref to process_unique_one
        if out1 and out2 == None:
            out = out1
        elif out2 and out1 == None:
            out = out2
        elif out1 and out2: # if r2 has hit both on top_ref and bottom_ref, run probability_hit
            choose = probability_hit2(out1[1], out2[1], distance)
            if choose == 0:
                out1[1][1] = "p"
                out = out1
            elif choose == 1:
                out2[1][1] = "p"
                out = out2
        else:
            out = [[get_output_var(hits[0][0]), 's'], []]

    elif len(hits[1]) == 1 and len(hits[2]) > 1 and len(hits[3]) > 1: # r1.CT against bottom_ref unique hit, r2 non-unique hit
        h1 = [hits[1], hits[2]]
        h2 = [hits[1], hits[3]]
        out1 = process_unique_one(h1, read1, read2, a, top_a, distance)
        out2 = process_unique_one(h2, read1, read2, a, bottom_a, distance)
        if out1 and out2 == None:
            out = out1
        elif out2 and out1 == None:
            out = out2
        elif out1 and out2:
            choose = probability_hit2(out1[1], out2[1], distance)
            if choose == 0:
                out1[1][1] = "p"
                out = out1
            elif choose == 1:
                out2[1][1] = "p"
                out = out2
        else:
            out = [[get_output_var(hits[1][0]), 's'], []]

    elif len(hits[2]) == 1 and len(hits[0]) > 1 and len(hits[1]) > 1: # r2.GA against top_ref unique hit, r1 non-unique hit
        h1 = [hits[0], hits[2]]
        h2 = [hits[1], hits[2]]
        out1 = process_unique_one(h1, read1, read2, a, top_a, distance)
        out2 = process_unique_one(h2, read1, read2, a, bottom_a, distance)
        if out1 and out2 == None:
            out = out1
        elif out2 and out1 == None:
            out = out2
        elif out1 and out2:
            choose = probability_hit2(out1[0], out2[0], distance)
            if choose == 0:
                out1[0][1] = "p"
                out = out1
            elif choose == 1:
                out2[0][1] = "p"
                out = out2
        else:
            out = [[], [get_output_var(hits[2][0]), 's']]

    elif len(hits[3]) == 1 and len(hits[0]) > 1 and len(hits[1]) > 1: # r2.GA against bottom_ref unique hit, r1 non-unique hit
        h1 = [hits[0], hits[3]]
        h2 = [hits[1], hits[3]]
        out1 = process_unique_one(h1, read1, read2, a, top_a, distance)
        out2 = process_unique_one(h2, read1, read2, a, bottom_a, distance)
        if out1 and out2 == None:
            out = out1
        elif out2 and out1 == None:
            out = out2
        elif out1 and out2:
            choose = probability_hit2(out1[0], out2[0], distance)
            if choose == 0:
                out1[0][1] = "p"
                out = out1
            elif choose == 1:
                out2[0][1] = "p"
                out = out2
        else:
            out = [[], [get_output_var(hits[3][0]), 's']]

    else:
        out = None

    return (out)


def process_unique_one(hits, read1, read2, a, converted_a, distance):
    """
    Method to process non unique mapping hits
    """
    if len(hits[0]) == 1 and len(hits[1]) > 0: # r1 is unique hit
        indx = 0
        seq = read2.GA # re-alignment sequence
        ref = converted_a # converted genome for re-alignment
    elif len(hits[1]) == 1 and len(hits[0]) > 0: # r2 is unique hit
        indx = 1
        seq = read1.CT
        ref = converted_a

    # Get reference sequence for unique hit
    local_ref = ref.seq(hits[indx][0].ctg, hits[indx][0].r_st - distance, hits[indx][0].r_en + distance)

    if local_ref: # build an index for local reference
        b = mp.Aligner(seq=local_ref, preset="sr", n_threads=1)
        newhits = []
        for hit in b.map(seq):
            if hit.mlen / (len(seq) * 1.0) > 0.8:
                newhits.append(hit)

        if len(newhits) == 1: # only 1 hit within local reference
            if indx == 0:
                out = [[get_output_var(hits[0][0]), "u"], [get_output_var(newhits[0]), "r"]]
                out[1][0]["ctg"] = out[0][0]["ctg"]
                out[1][0]["r_st"] = out[0][0]["r_st"] - distance + out[1][0]["r_st"]
                out[1][0]["r_en"] = out[0][0]["r_st"] - distance + out[1][0]["r_en"]
            else:
                out = [[get_output_var(newhits[0]), "r"], [get_output_var(hits[1][0]), "u"]]
                out[0][0]["ctg"] = out[1][0]["ctg"]
                out[0][0]["r_st"] = out[1][0]["r_st"] - distance + out[0][0]["r_st"]
                out[0][0]["r_en"] = out[1][0]["r_st"] - distance + out[0][0]["r_en"]

        elif len(newhits) > 0: # > 1 hits within local reference
            choose = probability_hit(newhits, distance) # According model, random choose
            if indx == 0:
                out = [[get_output_var(hits[0][0]), "u"], [get_output_var(newhits[choose]), "p"]]
                out[1][0]["ctg"] = out[0][0]["ctg"]
                out[1][0]["r_st"] = out[0][0]["r_st"] - distance + out[1][0]["r_st"]
                out[1][0]["r_en"] = out[0][0]["r_st"] - distance + out[1][0]["r_en"]
            else:
                out = [[get_output_var(newhits[choose]), "p"], [get_output_var(hits[1][0]), "u"]]
                out[0][0]["ctg"] = out[1][0]["ctg"]
                out[0][0]["r_st"] = out[1][0]["r_st"] - distance + out[0][0]["r_st"]
                out[0][0]["r_en"] = out[1][0]["r_st"] - distance + out[0][0]["r_en"]
        else:
            out = None

        return(out)


def probability_hit(newhits, distance): # random choose for > 1 hits within local reference
    distance_list = []
    for new_hit in newhits:
        if new_hit.r_st >= distance:
            distance_list.append(new_hit.r_st - distance)
        else:
            distance_list.append(distance - new_hit.r_en)
    cumulative_probability = 0
    probability_list = []
    for dist in distance_list:
        scale_probability = math.exp(-0.8 * dist / 50 - 0.6618) # the model
        if scale_probability < 0:
            scale_probability = 0
        cumulative_probability += scale_probability
        probability_list.append(cumulative_probability)

    selected_probability = random.random() * cumulative_probability # random choose hit according the probability of each hit

    choose = 0
    for a in range(len(probability_list) - 1): #
        if probability_list[a] <= selected_probability <= probability_list[a + 1]:
            choose = a + 1
    return (choose)


def probability_hit2(out1, out2, distance): # random choose for the read has hits both on top_ref and bottom_ref
    distance_list = []
    if out1[0]["r_st"] >= distance:
        distance_list.append(out1[0]["r_st"] - distance)
    else:
        distance_list.append(distance - out1[0]["r_en"])

    if out2[0]["r_st"] >= distance:
        distance_list.append(out2[0]["r_st"] - distance)
    else:
        distance_list.append(distance - out2[0]["r_en"])

    cumulative_probability = 0
    probability_list = []
    for dist in distance_list:
        scale_probability = math.exp(-0.8 * dist / 50 - 0.6618)
        if scale_probability < 0:
            scale_probability = 0
        cumulative_probability += scale_probability
        probability_list.append(cumulative_probability)

    selected_probability = random.random() * cumulative_probability

    choose = 0
    for a in range(len(probability_list) - 1):
        if probability_list[a] <= selected_probability <= probability_list[a + 1]:
            choose = a + 1
    return (choose)
