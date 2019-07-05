import multiproc
from sys import argv
import mappy as mp
import processReads
import alignment
import output

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

def printHelp():
    print("Usage: patcher -g <target.fa> -r1 <query.fq.gz> -r2 <query.fq.gz> -o <output.txt>")
    print("Options:")
    print("\t-g\t\tSet the reference sequence (Must be FASTA)")
    print("\t-r1\t\tSet the read sequence for 1 end (Must be FASTQ/Can be gzipped)")
    print("\t-r2\t\tSet the read sequence for 2 end (Must be FASTQ/Can be gzipped)")
    print("\t-t\t\tSet the number of threads")
    print("\t-d\t\tSet distance to map around total is +/- value")
    print("\t-o\t\tSet output filename")
    print("\t-c\t\tCut Site. Defaults to GATC. This will cut GATCGATC")
    print("\t-l\t\tMinumum length read to keep for mapping. Defaults to 20")
    exit()

    
def runsingle(ref,reads1,reads2,fname,distance,cutsite,minlen):
    fhandle = open(fname, "w")
    print("Loading in Reference")
    a = mp.Aligner(ref,preset="sr")
    if not a: raise Exception("ERROR: failed to load/build index file")
    print("Done.")
    print("Running Alignment")
    while True:
        try:
            r1=processReads.read(reads1.__next__())
            r1.splitRead(cutsite,minlen)
            r1.qualtrim(10,10)
            r2=processReads.read(reads2.__next__())
            r2.splitRead(cutsite,minlen)
            r2.qualtrim(10,10)
            #r1 and r2 could be None at this stage
            if r1.seq and r2.seq:
                res=alignment.map(r1,r2,a,1000)
                if res:
                    out=output.processOutput(res,r1,r2,"SAM",1000)
                    if out:
                        fhandle.write(f"{out[0]}\n")
                        fhandle.write(f"{out[1]}\n")                
        except StopIteration:
            break
    fhandle.close()


if __name__ == "__main__":
    distance=10
    nthreads=1
    cutsite="GATC"
    minlen=20
    myargs = getopts(argv)
    if '-g' in myargs: 
        ref=myargs["-g"]
    else:
        printHelp()
    if '-o' in myargs: 
        fname=myargs["-o"]
    else:
        printHelp()
    if '-r1' in myargs: 
        r1name=myargs["-r1"]
        reads1=mp.fastx_read(r1name)
    else:
        printHelp()
    if '-r2' in myargs: 
        r2name=myargs["-r2"]
        reads2=mp.fastx_read(r2name)
    else:
        printHelp()
    if '-d' in myargs:
        distance=int(myargs["-d"])
    if '-t' in myargs:
        nthreads=int(myargs["-t"])
    if '-c' in myargs:
        cutsite=myargs["-c"]
    if '-l' in myargs:
        minlen=int(myargs["-l"])
    
    
    print(f"Using refrence:{ref}")
    print(f"Using read 1:{r1name}")
    print(f'Using read 2:{r2name}')
    print(f'Using distance +/-:{distance}')
    print(f"Using threads:{nthreads}")
    print(f"Using cutsite:{cutsite}")
    print(f"Writing to:{fname}")
    if nthreads>3:
        multiproc.run(ref,reads1,reads2,fname,distance,nthreads,cutsite,minlen)
    else:
        if nthreads>1:
            print("Cannot run multithreading with less than 3 threads defaulting to single")
        runsingle(ref,reads1,reads2,fname,distance,cutsite,minlen)
    reads1.close()
    reads2.close()        
