import threading
import mappy as mp
import time
from multiprocessing import Process, Queue
from queue import Empty
import process_reads
import alignment
from output import SAMBAMWriter

def feed(queue, parlist):
    """
    Method to read data and feed into the calc queue
    """
    while True:
        try:
            if queue.full():
                time.sleep(2)
            else:
                read1 = process_reads.Read(parlist[0].__next__())
                read2 = process_reads.Read(parlist[1].__next__())
                queue.put((read1, read2))
        except StopIteration:
            print("Finished Feeding Data")
            for t in range(parlist[2]):
                queue.put(None)
            break

def process_read(queue_in, queue_out, thread_id, reference, distance, cut_site, min_len, debug):
    """
    Method for the main processing including quality trim and read splitting
    plus alignment
    """
    nreads = 0
    time_processing_reads = 0
    alignment_time = 0
    waiting_for_reads = 0
    while True:
        try:
            reads = queue_in.get(block=False)
            start = time.time()
            if reads:
                reads[0].split_read(cut_site, min_len)
                reads[0].qual_trim(10, 10)
                reads[1].split_read(cut_site, min_len)
                reads[1].qual_trim(10, 10)
                time_processing_reads += time.time() - start
                nreads += 1
                if nreads % 10000 == 0 and debug:
                    print(f"""
Time taken to process reads in: {time_processing_reads}
Time taken for alignment: {alignment_time}
Time taken waiting to process reads: {waiting_for_reads}
                          """)
                if reads[0].seq and reads[1].seq:
                    start = time.time()
                    res = alignment.map_reads(reference, reads[0], reads[1], distance)
                    alignment_time += time.time()-start
                    if res:
                        queue_out.put((res, reads[0], reads[1], distance))
            else:
                queue_out.put(None)
                break
        except Empty:
            start = time.time()
            time.sleep(0.1)
            waiting_for_reads += time.time()-start

def save(queue, fname, reference, output_type, nthreads):
    """
    Method to write out data
    """
    sambam_output = SAMBAMWriter(fname, reference, output_type)
    n_nones = 0
    while True:
        try:
            out = queue.get(block=False)
            if out:
                sambam_output.process_output(out[0], out[1], out[2], out[3])
            else:
                if out == None:
                    n_nones += 1
        except Empty:
            time.sleep(1)
            if n_nones == (nthreads - 2):
                break


def run(reffile, read1, read2, fname, distance, nthreads, cut_site, min_len, output_type, debug):
    """
    Method which sets up the multiprocessing and orchestrates the process
    """
    print("Loading in Reference")
    reference = mp.Aligner(reffile, preset="sr")
    if not reference:
        raise Exception("ERROR: failed to load/build index file")
    print("Done.")

    read_queue = Queue(50000)
    result_queue = Queue()
    parlist = (read1, read2, nthreads)
    feed_proc = Process(target=feed, args=(read_queue, parlist))
    writ_proc = Process(target=save, args=(result_queue, fname, reference, output_type, nthreads))

    feed_proc.start()

    threads = []
    for index in range(nthreads - 2):
        x = Process(target=process_read,
                             args=(read_queue, result_queue, index, reference, distance, cut_site, min_len, debug),
                             daemon=True)
        x.start()
        threads.append(x)

    writ_proc.start()

    feed_proc.join()
    for index, thread in enumerate(threads):
        thread.join()
    writ_proc.join()
