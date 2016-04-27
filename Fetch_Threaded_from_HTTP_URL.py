#!/usr/bin/env python
import Queue
import threading
import urllib2
import time,sys

hosts =[] #["http://yahoo.com", "http://google.com", "http://amazon.com","http://ibm.com", "http://apple.com"]
for each in open("/media/srikanth/My Passport/CellMap_RNASeq/files_cellmap"):
	ident=re.search("JH-(\d+).*L(\d+)_R[12]_001.fastq.gz",each.strip()).group(1)
	if ident==i:
		hosts.append(each.strip())
		
queue = Queue.Queue()

class ThreadUrl(threading.Thread):
    """Threaded Url Grab"""
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            #grabs host from queue
            host = self.queue.get()

            #grabs urls of hosts and prints first 1024 bytes of page
            
			urllib.urlretrieve("http://genome.nyumc.org/results/external/JHU/Burns_HSQ_2/2014-08-19/fastq/"+host,host )
            print "Downloading", host


            #signals to queue job is done
            self.queue.task_done()

start = time.time()
def main():
    
    #spawn a pool of threads, and pass them queue instance 
    for i in range(6):
        t = ThreadUrl(queue)
        t.setDaemon(True)
        t.start()

    #populate queue with data
    for host in hosts:
        queue.put(host)
    
    #wait on the queue until everything has been processed
    queue.join()
main()
print "Elapsed Time: %s" % (time.time() - start)
