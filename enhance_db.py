"""
This script was created to enhance the data gathered in Inovatoxin.
It extracts more relevant data from SwissProt and updates the DB.
"""

import datetime
import threading
# Possible errors
from urllib.error import HTTPError, URLError

# Retrieve Protein Names
from Bio import ExPASy, SeqIO
# Nice Warnings
from colorama import Back, Fore, Style
# Connect to MongoDB
from pymongo import MongoClient

from connect import connectInovatoxin

# Lock to sync threads
LOCK = threading.Lock()
BAR = threading.Lock()
SEM = threading.Semaphore(value=4)

# Number of threads
WORKERS = 12

class GeneUpdater(threading.Thread):
    def __init__(self, work_list, id):
        threading.Thread.__init__(self)
        self.genes = work_list
        self.id = id
    def run(self):
        global PROGRESS
        my_data = []
        for g in self.genes:
            uni_number = get_uniprot_number(g)
            BAR.acquire()
            PROGRESS.update(g)
            BAR.release()
            if uni_number is not None:
                my_data.append({'gene':g, 'id':uni_number})
        db = connectInovatoxin()
        print(Fore.BLUE + "Worker %d moving to phase 2" % self.id + Style.RESET_ALL)
        for p in my_data:
            update_entry(db, p)
        print("\nWorker %d exiting" % self.id)

class ProgressBar():
    """ Used to display progress when collecting genes.
    """
    def __init__(self):
        """ Instantiate progress bar"""
        self.total = 0
        self.start_time = 0
        self.done = 0
        print("ProgressBar created")

    def start(self, total_genes):
        """ Starts progress bar """
        self.total = total_genes
        self.start_time = datetime.datetime.now()
        self.done = 0
        print("CURRENT GENE -  TIME ELAPSED  -  AVERAGE TIME  - CURRENT GENE")
  
    def update(self, gene):
        """ Updates screen progress bar"""
        self.done += 1
        now = datetime.datetime.now()
        avg = (now - self.start_time)/self.done
        print("\r [{:04}/{:04}] - {} - {} - Current protein: {:12s} ".format(self.done, self.total,
            now - self.start_time, avg, gene), end="")

def get_uniprot_number(gene, retry=0):
    """ Returns information from SwissProt on gene """
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            info = SeqIO.read(handle, "swiss")
            return info.id
    except:
        if(retry < 10):
            print("\nConnection failed for gene %s. Retrying... (%d)" % (gene, retry))
            return get_uniprot_number(gene, retry+1)
    print(Fore.YELLOW + "WARNING: " + Style.RESET_ALL + "%s not found" % gene)
    return None

def update_entry(db, data):
    """ Updates an entry in MongoDB """
    SEM.acquire()
    db.update_many({'Blast_fullAccession':data['gene']}, {'$set': {
        'Uniprot_accession': data['id'],
    }})
    SEM.release()

PROGRESS = ProgressBar()

def main():
    """ Aggregates Mongo by gene name. Updates values """
    global PROGRESS
    db = connectInovatoxin()
    genes = db.distinct('Uniprot_accession')
    size = len(genes)
    # Start progress bar
    PROGRESS.start(size)
    # Prepare workers duties
    worker_part = size // WORKERS
    threads = []
    for i in range(WORKERS):
        sts = i * worker_part
        end = sts + worker_part
        if i < WORKERS-1:
            t = GeneUpdater(genes[sts : end], i)
        else:
            t = GeneUpdater(genes[sts:], i)
        t.start()
        threads.append(t)
    # Join with workers
    for t in threads:
        t.join()
    # Finishes
    del(PROGRESS)
    print( Fore.GREEN + "All done!" + Style.RESET_ALL)
        

if __name__ == '__main__':
    main()
