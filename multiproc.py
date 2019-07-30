from queue import Empty
import processReads
import alignment
import output
import mappy as mp
from multiprocessing import Process, Queue

def feed(queue, parlist):
    while True:
        try:
            r1=processReads.read(parlist[0].__next__())
            r2=processReads.read(parlist[1].__next__())
            queue.put((r1,r2))
        except StopIteration:
            print("Finished Feeding Data")
            break

def calc(queueIn, queueOut,threadid,d,cutsite,minlen):
    nreads=0
    nNone=0
    while True:
        try:
            r1,r2 = queueIn.get(block = False)
            if r1 and r2:
                nNone=0
                r1.splitRead(cutsite,minlen)
                r1.qualtrim(10,10)
                r2.splitRead(cutsite,minlen)
                r2.qualtrim(10,10)
                nreads+=1
                if r1.seq and r2.seq:
                    res=alignment.map(r1,r2,a,d)
                    if res:
                        queueOut.put(res,r1,r2,d)
        except Empty:
            nNone+=1
            if nNone>1000000:
                queueOut.put(None)
                print(f"Finshed Working Thread:{threadid}")
                break

def write(queue, fname,nthreads):


    fhandle = open(fname, "w")
    nNones=0
    while True:
        try:
            out = queue.get(block = False)
            if out:
                out = output.processOutput(out[0],out[1],out[2],out[3])               
            else:
                if out==None:
                    nNones+=1
        except Empty:
            if nNones==(nthreads-2):
                break
    fhandle.close()

def run(ref,r1,r2,sambam_output,d,nthreads,cutsite,minlen):
    global a
    print("Loading in Reference")
    a = mp.Aligner(ref,preset="sr")
    if not a: raise Exception("ERROR: failed to load/build index file")
    print("Done.")
    print("Running Alignment")
    workerQueue = Queue()
    writerQueue = Queue()
    parlist = (r1,r2)
    feedProc = Process(target = feed , args = (workerQueue, parlist))
    calcProc = [Process(target = calc , args = (workerQueue, writerQueue,i,d,cutsite,minlen)) for i in range(nthreads-2)]
    writProc = Process(target = write, args = (writerQueue, sambam_output,nthreads))


    feedProc.start()
    for p in calcProc:
        p.start()
    writProc.start()

    feedProc.join ()
    for p in calcProc:
        p.join()
    writProc.join ()
    print(f"Finished")
