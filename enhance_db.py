"""
This script was created to enhance the data gathered in Inovatoxin.
It extracts more relevant data from SwissProt and updates the DB.
"""
# Possible errors
from urllib.error import HTTPError

# Retrieve Protein Names
from Bio import ExPASy, SeqIO
# Connect to MongoDB
from pymongo import MongoClient
# Nice Warnings
from colorama import Fore, Back, Style


def connectInovatoxin():
    """ Connection to MongoDB """
    client = MongoClient()
    db = client['inovatoxin']
    return db['Proteins_protein']

def getSwissProtInfo(gene):
    """ Returns information from SwissProt on gene """
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            return SeqIO.read(handle, "swiss")
    except HTTPError:
        print(Fore.YELLOW + "WARNING: " + Style.RESET_ALL + "%s not found" % gene)
    return None

def getExtraInfoForEntry(gene):
    """ Updates the entry adding extra information """
    info = getSwissProtInfo(gene)
    data = {}
    data['id'] = info.id
    data['gene'] = gene
    for db in info.dbxrefs:
        if db.startswith('KEGG:'):
            data['kegg'] = db[5:]
        elif db.startswith('STRING:'):
            data['string'] = db[7:] 
        elif db.startswith('eggNOG:'):
            data['eggnog'] = db[7:]
        elif db.startswith('KO:'):
            data['ko'] = db[3:]
    return data

def UpdateEntry(data):
    """ Updates an entry in MongoDB """
    db = connectInovatoxin()
    db.update_many({'Blast_fullAccession':data['gene']}, {'$set': {
        'Uniprot_accession': data['id'],
        'KEGG_org': data.get('kegg', None),
        'KEGG_ko': data.get('ko', None),
        'STRING': data.get('string', None),
        'EggNOG_indexTerm': data.get('eggnog', None)
        }}
    )

def main():
    """ Aggregates Mongo by gene name. Updates values """
    db = connectInovatoxin()
    genes = db.distinct('Uniprot_accession')
    size = len(genes)
    done = 0
    for gene in genes:
        print('[{}/{}] - {}\r'.format(done, size, gene), end='')
        info = getExtraInfoForEntry(gene)
        if info is not None:
            UpdateEntry(info)
        done += 1

if __name__ == '__main__':
    main()
