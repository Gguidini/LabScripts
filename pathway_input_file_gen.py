""" This script was created to search for EC reference in Uniprot,
    gather this info with others in Inovatoxin database and generate files.

    THe files generated are the input files for Pathway Tools.
"""
import re
import os
# Possible errors
from urllib.error import HTTPError

# Retrieve Protein Names
from Bio import ExPASy, SeqIO
# Connect to MongoDB
from pymongo import MongoClient
# Nice Warnings
from colorama import Fore, Back, Style

def get_protein_EC(gene):
    rgx = re.compile(r"EC=\d+\.\d+\.\d+\.\d+")
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            seq_record = SeqIO.read(handle, "swiss")
            match = rgx.search(seq_record.description)
            if match is not None:
                return match.group(0)
    except HTTPError:
        print(Fore.YELLOW + "WARNING: " + Style.RESET_ALL + "%s not found" % gene)
    return None

def connect_Inovatoxin():
    """ Connection to MongoDB """
    client = MongoClient()
    MONGO = client['inovatoxin']
    return MONGO['Proteins_protein']


DB = connect_Inovatoxin()
DOCS = DB.find()
COUNT = DB.count_documents({})
i = 1

DIR = os.getcwd()
try:
    os.mkdir("fastas")
except FileExistsError:
    print("Dir fastas exists. Skipping")
try:
    os.mkdir("infos")
except FileExistsError:
    print("Dir infos exists. Skipping")

GENETIC_ELEMENTS = open("genetic_elements.dat", 'w')

GENES = DB.distinct("Blast_fullAccession")
gene_map = {}
for g in GENES:
    print("\r[{}/{}] - Getting EC for {}".format(i, len(GENES), g), end="")
    i += 1
    gene_map[g] = get_protein_EC(g)

i = 1
for doc in DOCS:
    print("\r[{}/{}] - {}".format(i, COUNT, doc["Blast_fullAccession"]), end="")
    i += 1

    ec = gene_map[doc["Blast_fullAccession"]]
    if ec is None:
        print(Fore.RED + " EC NOT FOUND" + Style.RESET_ALL)
        continue

    prot_id = "Inova_" + doc["Blast_fullAccession"] + "_" + str(doc["id"])
    fasta_file = open("fastas/" + prot_id + ".fsa", 'w')
    info_file = open("infos/" + prot_id + ".pf", 'w')
    seq = doc["transcript_sequence"]

    # General file
    text = "ID  " + prot_id + "\n"
    text += "NAME  " + doc["name"] + "\n"
    text += "TYPE  :CONTIG\n"
    text += "CIRCULAR?  N\n"
    text += "SEQ-FILE    " + DIR + "fastas/" + prot_id + ".fsa\n"
    text += "ANNOT-FILE    " + DIR + "infos/" + prot_id + ".pf\n"
    text += "//\n"
    GENETIC_ELEMENTS.write(text)

    # Fasta file
    fasta_file.write(">" + prot_id + "\n")
    fasta_file.write(seq)
    fasta_file.close()

    # info file
    text = "ID  " + prot_id + "\n"
    text += "NAME  " + doc["name"] + "\n"
    text += "PRODUCT-TYPE  P\n"
    text += "STARTBASE  1\n"
    text += "ENDBASE  " + str(len(seq)) + "\n"
    text += "EC  " + ec + "\n"
    text += "//\n"
    info_file.write(text)
    info_file.close()

