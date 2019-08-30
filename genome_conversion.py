"""
Module containing the function for genome preparation
"""
def conversion(ref,base1,base2):
    new_file = "Reference" + "_" + base1 + "_" + base2 + ".fasta"
    fle=open(new_file,"w")
    for i in ref:
        name = i[0]
        seq = i[1].replace(base1, base2)
        seq = seq.replace(base1.lower(), base2.lower())
        fle.write(f">{name}\n")
        fle.write(f"{seq}\n")
    fle.close()