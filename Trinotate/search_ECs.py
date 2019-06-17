# This script queries EC numbers for proteins based on Uniprot name

import datetime
import sqlite3
import sys
import threading
import re
import time
from urllib.error import HTTPError, URLError

import colorama
import requests
from Bio import ExPASy, SeqIO
from bs4 import BeautifulSoup

########################################################### CONSTANTS AND SHARED VARIABLES
NUMBER_OF_WORKERS = 20
DATABASE = ""
if len(sys.argv) == 1:
    DATABASE = input("Path to SQLite file: ")
else:
    DATABASE = sys.argv[1]
CONN = sqlite3.connect(DATABASE)
CURSOR = CONN.cursor()

EC_QUERIES = []
EC_QUERIES_LOCK = threading.Lock()
WRITERS_LOCK = threading.Lock()

########################################################### CLASSES
class ProgressBar():
    """ Used to display progress when collecting genes.
    """
    def __init__(self):
        print("ProgressBar created")

    def start(self, total_genes):
        self.total = total_genes
        self.start_time = datetime.datetime.now()
        self.done = 0
        print("CURRENT GENE -  TIME ELAPSED  -  AVERAGE TIME  - CURRENT GENE")
   
    def update(self, gene):
        self.done += 1
        now = datetime.datetime.now()
        avg = (now - self.start_time)/self.done
        print("\r [{:04}/{:04}] - {} - {} - Getting EC for {:12s} ".format(self.done, self.total,
            now - self.start_time, avg, gene), end="")

class WorkerGetEC(threading.Thread):
    """ Speedup the work of searching EC references
    """
    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        my_conn = sqlite3.connect(DATABASE)
        my_cursor = my_conn.cursor()
        row = ""
        while row is not None:
            row = get_next_query()
            idx = get_protein_EC(row[0])
            WRITERS_LOCK.acquire()
            PROGRESS.update(row[0])
            if idx is not None:
                my_cursor.execute("INSERT INTO UniprotIndex('Accession', 'LinkId', 'AttributeType') VALUES (?,?,F);", (row[0], idx))
            WRITERS_LOCK.release()

########################################################### FUNCTIONS
def get_protein_EC(gene, retry=0):
    """ Queries Uniprot for a gene entry and extracts the EC, if any.
        If the gene is successfully queried, but no EC is present, returns None.
        It's possible that, due to connection problems, a gene that is in
        Uniprot is not found, so it will try again after a cooldown period.

        > Input
        gene : str => the gene code to be queried
        retry : int => number of tries. Max 10.

        > Output
        - EC for GENE, if GENE has one annotated in Uniprot.
        - None, if GENE doesn't have an EC
        - Exception, if any exception occurred.
          Most common exceptions are HTTPError or ValueError.
    """
    rgx = re.compile(r"EC=\d+\.\d+\.\d+\.\d+")
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            seq_record = SeqIO.read(handle, "swiss")
            match = rgx.search(seq_record.description)
            if match is not None:
                return match.group(0)
    except Exception as e:
        if retry < 10:
            time.sleep(5) # cool down time 5s
            print("\nGENE NOT FOUND. RETRYING (%d)" % retry)
            return get_protein_EC(gene, retry+1)
        return e
    except KeyboardInterrupt as k:
        print("\nKeyBoard Interrupt Signal received. Aborting")
        return k
    return None

def get_next_query():
    """ Gets the next Uniprot index to be searched for AS reference.
        Returns None if there's no more work
    """
    EC_QUERIES_LOCK.acquire()
    if EC_QUERIES != []:
        ret = EC_QUERIES.pop(0)
    else:
        ret = None
    EC_QUERIES_LOCK.release()
    return ret

PROGRESS = ProgressBar()
if __name__ == "__main__":
    CURSOR.execute('SELECT DISTINCT(BlastDbase.FullAccession) FROM BlastDbase;')
    EC_QUERIES = CURSOR.fetchall()
    PROGRESS.start(len(EC_QUERIES))
    workers = []
    for i in range(NUMBER_OF_WORKERS):
        w = WorkerGetEC()
        w.start()
        workers.append(w)
    # Wait for work to be completed
    for w in workers:
        w.join()
    CONN.commit()
    CONN.close()
    