# This script cleans SQLite database data.
# The databases where created using Trinotate
# Blasts against cusom databases didn't go so well into the tables.

import datetime
import sqlite3
import sys
import threading
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

ARACH_INDEX_QUERIES = []
ARACH_QUERY_LOCK = threading.Lock()





########################################################## CLASSES
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

class WorkerUniprot2Arach(threading.Thread):
    """ Speedup the work of searching Uniprot for Arachnoserver references
    """
    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        row = ""
        while row is not None:
            row = get_next_query_arach()
            PROGRESS.update(row[1])
            idx = get_arachindex_from_uniprot(row[1])
            if idx is not None:
                CURSOR.execute("UPDATE BlastDbase SET ArachnoserverIndex = ? WHERE TrinityID = ? AND DatabaseSource != 'arachnoserver.pep.fa';", (idx, row[0]))

def get_next_query_arach():
    """ Gets the next Uniprot index to be searched for AS reference.
        Returns None if there's no more work
    """
    ARACH_QUERY_LOCK.acquire()
    if ARACH_INDEX_QUERIES != []:
        ret = ARACH_INDEX_QUERIES.pop(0)
    else:
        ret = None
    ARACH_QUERY_LOCK.release()
    return ret

def get_arachnoserver_name(idx):
    """ Gets the Arachnoserver protein name from the index number"""
    page = requests.get("http://www.arachnoserver.org/toxincard.html", params={'id':idx})
    soup = BeautifulSoup(page.text, 'html.parser')
    name = soup.find(class_='toxinname')
    return name.text

def get_arachnoserver_index(name):
    """ Gets the Arachnoserver index number from the name"""
    page = requests.get("http://www.arachnoserver.org/basicsearch.html", params={'keywords':name})
    soup = BeautifulSoup(page.text, 'html.parser')
    carddata = soup.find(class_='carddata')
    links = carddata.find_all('a', href=True)
    n = links[0]['href']
    return n[n.find('?id=')+4:]

def get_arachindex_from_uniprot(gene, retry=0):
    """ Gets Arachnoserver index from Uniprot, if any"""
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            info = SeqIO.read(handle, "swiss")
            for db in info.dbxrefs:
                if db.startswith('ArachnoServer'):
                    return db[14:]
            return None
    except:
        if(retry < 10):
            print("\nConnection failed for gene %s. Retrying... (%d)" % (gene, retry))
            return get_uniprot_number(gene, retry+1)
    print(colorama.Fore.YELLOW + "WARNING: " + colorama.Style.RESET_ALL + "%s not found" % gene)
    return None

def get_uniprot_number(gene, retry=0):
    """ Gets the name of Uniprot entry from index"""
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            info = SeqIO.read(handle, "swiss")
            return info.name
    except HTTPError:
        return None
    except URLError:
        if(retry < 10):
            print("\nConnection failed for gene %s. Retrying... (%d)" % (gene, retry))
            return get_uniprot_number(gene, retry+1)
    print(colorama.Fore.YELLOW + "WARNING: " + colorama.Style.RESET_ALL + "%s not found" % gene)
    return None

def clean_warns():
    """Remove ROWS that have script output in them.
        These rows are not usefull in any way.
    """
    print("Checking for Rows with Warning as TrinityID....", end="")
    CURSOR.execute('SELECT count(*) FROM BlastDbase WHERE BlastDbase.TrinityID = "Warning:" GROUP BY BlastDbase.TrinityID;')
    if CURSOR.fetchone() is not None:
        print("yes")
        CURSOR.execute('DELETE FROM BlastDbase WHERE BlastDbase.TrinityID = "Warning:"')
        print("%d rows removed by Warning clean-up." % CURSOR.rowcount)
    else:
        print("no")

def clean_toxprot_accession():
    """ Cleans ToxProt accession format"""
    all_toxprot_entries = """
    SELECT
        DISTINCT(BlastDbase.TrinityID),
        BlastDbase.UniprotSearchString,
        BlastDbase.FullAccession,
        BlastDbase.DatabaseSource
    FROM
        BlastDbase
    WHERE
        BlastDbase.DatabaseSource = "toxprot";
    """
    toxprot_to_query = open("toxprot_entries_to_query.txt", "w")
    for row in CONN.cursor().execute(all_toxprot_entries):
        row_info = row[1].split('|')
        if len(row_info) == 3:
            CURSOR.execute('SELECT count(*) FROM UniprotIndex WHERE UniprotIndex.Accession = ? GROUP BY UniprotIndex.Accession;', (row_info[2],))
        if len(row_info) == 3 and CURSOR.fetchone()[0] is None:
            print("%s (%s) not found in UniprotIndex" % (str(row[0]), str(row_info[2])))
            toxprot_to_query.write(row[0], row_info[2], "\n")
        elif len(row_info) == 3:
            CURSOR.execute('UPDATE BlastDbase SET FullAccession = ?, UniprotSearchString = ? WHERE TrinityID = ? AND DatabaseSource = "toxprot";', (row_info[2], row_info[2], row[0]))
        else:
            print(colorama.Fore.YELLOW + "Strange ROW:" + colorama.Style.RESET_ALL, row)
    toxprot_to_query.close()

def has_arachnoserver_column():
    """ Verifies if ArachnoserverIndex columns has been created """
    for r in CURSOR.execute('pragma table_info(BlastDbase)'):
        if r[1] == 'ArachnoserverIndex':
            return True
    return False

def has_arachnoserver_table():
    """ Verifies if ArachnoserverReference table has been created """
    CURSOR.execute('SELECT name FROM sqlite_master WHERE type="table" AND name="ArachnoserverReference";')
    if CURSOR.fetchone() is not None:
        return True
    return False

def clean_arachnoserver():
    """ Cleans up Arachnoserver index entries"""
    if not has_arachnoserver_column():
        CURSOR.execute("ALTER TABLE BlastDbase ADD COLUMN ArachnoserverIndex TEXT;")
    if not has_arachnoserver_table():
        create_table = """
            CREATE TABLE "ArachnoserverReference" (
            "Index"	INTEGER,
            "Name"	TEXT,
            PRIMARY KEY("Index")
        );
        """
        CURSOR.execute(create_table)
    ALL_ARACHNOSERVER = """
        SELECT
            DISTINCT(BlastDbase.TrinityID),
            BlastDbase.UniprotSearchString,
            BlastDbase.FullAccession,
            BlastDbase.DatabaseSource
        FROM
            BlastDbase
        WHERE
            BlastDbase.DatabaseSource = "arachnoserver.pep.fa";
    """
    for row in CONN.cursor().execute(ALL_ARACHNOSERVER):
        row_info = row[1].split('|')
        if len(row_info) == 3:
            row_info[0] = row_info[0].replace('as:', "")
            if(row_info[1].startswith('sp:')):
                row_info[1] = row_info[1].replace('sp:', '')
                name = get_uniprot_number(row_info[1])
                CURSOR.execute("UPDATE BlastDbase SET FullAccession = ?, UniprotSearchString = ?, ArachnoserverIndex = ? WHERE TrinityID = ? AND DatabaseSource = 'arachnoserver.pep.fa';",
                (name, name, row_info[2], row[0]))
                try:
                    CURSOR.execute("INSERT INTO ArachnoserverReference('Index', 'Name') VALUES (?,?);", (int(row_info[2]), row_info[0]))
                except sqlite3.IntegrityError:
                    pass # Entry already there
            else:
                CURSOR.execute("UPDATE BlastDbase SET FullAccession = ?, UniprotSearchString = ?, ArachnoserverIndex = ? WHERE TrinityID = ? AND DatabaseSource = 'arachnoserver.pep.fa';",
                (None, None, row_info[2], row[0]))
                try:
                    CURSOR.execute("INSERT INTO ArachnoserverReference('Index', 'Name') VALUES (?,?);", (int(row_info[2]), row_info[0]))
                except sqlite3.IntegrityError:
                    pass
        elif len(row_info) == 2:
            row_info[0] = row_info[0].replace("as:", "")
            name = get_uniprot_number(row_info[0])
            if name is not None:
                CURSOR.execute("UPDATE BlastDbase SET FullAccession = ?, UniprotSearchString = ?, ArachnoserverIndex = ? WHERE TrinityID = ? AND DatabaseSource = 'arachnoserver.pep.fa';",
                (name, name, row_info[1], row[0]))
                try:
                    CURSOR.execute("INSERT INTO ArachnoserverReference('Index', 'Name') VALUES (?,?);", (int(row_info[1]), get_arachnoserver_name(row_info[1])))
                except sqlite3.IntegrityError:
                    pass # Entry already there
            else:
                CURSOR.execute("UPDATE BlastDbase SET FullAccession = ?, UniprotSearchString = ?, ArachnoserverIndex = ? WHERE TrinityID = ? AND DatabaseSource = 'arachnoserver.pep.fa';",
                (None, None, row_info[1], row[0]))
                try:
                    CURSOR.execute("INSERT INTO ArachnoserverReference('Index', 'Name') VALUES (?,?);", (int(row_info[1]), row_info[0]))
                except sqlite3.IntegrityError:
                    pass # Entry already there
        else:
            row_info[0] = row_info[0].replace('as:', '')
            idx = get_arachnoserver_index(row_info[0])
            CURSOR.execute("UPDATE BlastDbase SET FullAccession = ?, UniprotSearchString = ?, ArachnoserverIndex = ? WHERE TrinityID = ? AND DatabaseSource = 'arachnoserver.pep.fa';",
                (None, None, idx, row[0]))
            try:
                CURSOR.execute("INSERT INTO ArachnoserverReference('Index', 'Name') VALUES (?,?);", (int(idx), row_info[0]))
            except sqlite3.IntegrityError:
                pass # Entry already there

######################################################################### MAIN
PROGRESS = ProgressBar()
if __name__ == "__main__":
    try:
        clean_warns()
        CONN.commit()
        print(colorama.Fore.GREEN, "Cleaning Warnings successful.", colorama.Style.RESET_ALL)
    except sqlite3.OperationalError as exp:
        CONN.rollback()
        print(colorama.Fore.RED, "Cleaning Warnings created an error. Rolled back.", colorama.Style.RESET_ALL)
        print(exp)
    try:
        clean_toxprot_accession()
        CONN.commit()
        print(colorama.Fore.GREEN, "Cleaning toxprot successful.", colorama.Style.RESET_ALL)
    except sqlite3.OperationalError as exp:
        CONN.rollback()
        print(colorama.Fore.RED, "Cleaning toxprot created an error. Rolled back.", colorama.Style.RESET_ALL)
        print(exp)
    try:
        clean_arachnoserver()
        CONN.commit()
        print(colorama.Fore.GREEN, "Cleaning arachnoserver successful.", colorama.Style.RESET_ALL)
    except sqlite3.OperationalError as exp:
        CONN.rollback()
        print(colorama.Fore.RED, "Cleaning arachnoserver created an error. Rolled back.", colorama.Style.RESET_ALL)
        print(exp)
    all_not_arch = """
        SELECT
            distinct(BlastDbase.TrinityID),
            BlastDbase.UniprotSearchString
        FROM
            BlastDbase
        WHERE
            BlastDbase.DatabaseSource = 'toxprot' OR
            BlastDbase.DatabaseSource = 'Swissprot';
    """
    CURSOR.execute(all_not_arch)
    ARACH_INDEX_QUERIES = CURSOR.fetchall()
    PROGRESS.start(len(ARACH_INDEX_QUERIES))
    workers = []
    for i in range(NUMBER_OF_WORKERS):
        w = WorkerUniprot2Arach()
        w.start()
        workers.append(w)
    # Wait for work to be completed
    for w in workers:
        w.join()
    CONN.commit()
    CONN.close()   
