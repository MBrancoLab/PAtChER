def processOutput(res,read1,read2,outtype,distance):
    #TODO handle unpaired reads
    tag=[65,129]
    hit1=res[0][0]
    hit2=res[1][0]
    if outtype=="SAM":  
        tlen = "*"
        if res[0][1]=="u" and res[1][1]=="u":
            if hit1.ctg == hit2.ctg:
                if hit1.r_st<hit2.r_st:
                    tlen=hit2.r_en-hit1.r_st
                else:
                    tlen=hit1.r_en-hit2.r_st
            s1=f"{read1.name}\t{tag[0]}\t{hit1.ctg}\t{hit1.r_st}\t{hit1.mapq}\t{hit1.cigar_str}\t{hit2.ctg}\t{hit2.r_st}\t{tlen}\t{read1.seq}\t{read1.qual_str}\tPT:A:{res[0][1]}"
            s2=f"{read2.name}\t{tag[1]}\t{hit2.ctg}\t{hit2.r_st}\t{hit2.mapq}\t{hit2.cigar_str}\t{hit1.ctg}\t{hit1.r_st}\t{tlen}\t{read2.seq}\t{read2.qual_str}\tPT:A:{res[1][1]}"
        elif res[0][1]!="u":
            if hit1.r_en>distance:
                tlen = hit1.r_en - distance
            else:
                tlen = distance - hit1.r_st
            s1=f"{read1.name}\t{tag[0]}\t{hit2.ctg}\t{hit2.r_st-distance+hit1.r_st}\t{hit1.mapq}\t{hit1.cigar_str}\t{hit2.ctg}\t{hit2.r_st}\t{tlen}\t{read1.seq}\t{read1.qual_str}\tPT:A:{res[0][1]}"
            s2=f"{read2.name}\t{tag[1]}\t{hit2.ctg}\t{hit2.r_st}\t{hit2.mapq}\t{hit2.cigar_str}\t{hit2.ctg}\t{hit2.r_st-distance+hit1.r_st}\t{tlen}\t{read2.seq}\t{read2.qual_str}\tPT:A:{res[1][1]}"
        elif res[1][1]!="u":
            if hit2.r_en>distance:
                tlen = hit2.r_en - distance
            else:
                tlen = distance - hit2.r_st
            s1=f"{read1.name}\t{tag[0]}\t{hit1.ctg}\t{hit1.r_st}\t{hit1.mapq}\t{hit1.cigar_str}\t{hit1.ctg}\t{hit1.r_st-distance+hit2.r_st}\t{tlen}\t{read1.seq}\t{read1.qual_str}\tPT:A:{res[0][1]}"
            s2=f"{read2.name}\t{tag[1]}\t{hit1.ctg}\t{hit1.r_st-distance+hit2.r_st}\t{hit2.mapq}\t{hit2.cigar_str}\t{hit1.ctg}\t{hit1.r_st}\t{tlen}\t{read2.seq}\t{read2.qual_str}\tPT:A:{res[1][1]}"
        return((s1,s2))
