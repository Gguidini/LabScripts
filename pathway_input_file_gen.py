""" This script was created to search for EC reference in Uniprot,
    gather this info with others in Inovatoxin database and generate files.

    THe files generated are the input files for Pathway Tools.
"""
import re
import os
import datetime
import time
import pickle
# Possible errors
from urllib.error import HTTPError

# Retrieve Protein Names
from Bio import ExPASy, SeqIO
# Connect to MongoDB
from pymongo import MongoClient
# Nice Warnings
from colorama import Fore, Back, Style

def get_protein_EC(gene, retry=0):
    rgx = re.compile(r"EC=\d+\.\d+\.\d+\.\d+")
    try:
        with ExPASy.get_sprot_raw(gene) as handle:
            seq_record = SeqIO.read(handle, "swiss")
            match = rgx.search(seq_record.description)
            if match is not None:
                return match.group(0)
    except Exception as e:
        if retry < 10:
            time.sleep(.5) # cool down time
            return get_protein_EC(gene, retry+1)
        print(Fore.YELLOW + "WARNING: " + Style.RESET_ALL + "%s not found" % gene)
        return e
    return None

def connect_Inovatoxin():
    """ Connection to MongoDB """
    client = MongoClient()
    MONGO = client['inovatoxin']
    return MONGO['Proteins_protein']

def gen_files(gene_map, docs, species):
    # Information for progress feedback
    count = docs.count()
    it = 1
    # Current directory
    # Creating necessary directories

    try:
        os.mkdir(species + "/fastas")
    except FileExistsError:
        print("Dir fastas exists. Skipping")
    try:
        os.mkdir(species + "/infos")
    except FileExistsError:
        print("Dir infos exists. Skipping")
    # General file for species
    genetic_elements = open(species + '/genetic-elements.dat', 'w')
    # Generate files
    for doc in docs:
        print("\r[{}/{}] - {:12s}".format(it, count, doc["Blast_fullAccession"]), end="")
        it += 1

        ec = gene_map.get(doc["Blast_fullAccession"], None)
        if ec is None:
            print(Fore.RED + " EC NOT FOUND" + Style.RESET_ALL)
            continue
        ec = ec[3:]

        prot_id = "Inova_" + doc["Blast_fullAccession"] + "_" + str(doc["id"])
        fasta_file = open(species + "/fastas/" + prot_id + ".fsa", 'w')
        info_file = open(species + "/infos/" + prot_id + ".pf", 'w')
        seq = doc["transcript_sequence"]

        # General file
        text = "ID  " + prot_id + "\n"
        text += "NAME  " + doc["name"] + "\n"
        text += "TYPE  :CONTIG\n"
        text += "CIRCULAR?  N\n"
        text += "SEQ-FILE    " + species + "/fastas/" + prot_id + ".fsa\n"
        text += "ANNOT-FILE    " + species + "/infos/" + prot_id + ".pf\n"
        text += "//\n"
        genetic_elements.write(text)

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
        text += "DBLINK  " + doc['GO_id'] + '\n'
        text += "//\n"
        info_file.write(text)
        info_file.close()
    genetic_elements.close()

DB = connect_Inovatoxin()
i = 1
initial_time = datetime.datetime.now()
GENES = DB.distinct("Blast_fullAccession")
GENE_MAP = {}
print("CURRENT GENE - TIME ELAPSED - AVERAGE TIME - CURRENT GENE")
for g in GENES:
    time_before = datetime.datetime.now()
    print("\r[{}/{}] - {} - {} - Getting EC for {:12s} ".format(i, len(GENES),
        time_before - initial_time, ((time_before - initial_time)/i), g), end="")
    i += 1
    GENE_MAP[g] = get_protein_EC(g)
print(Fore.GREEN + "All ECs collected!" + Style.RESET_ALL)
# Gets docs separated by species 
SPIDER = DB.find({"has_spider": {"$gt": 0}})
try: 
    os.mkdir('spider')
except:
    print("spider exists, skipping")
SCORPION = DB.find({"has_scorpion": {"$gt": 0}})
try: 
    os.mkdir('scorpion')
except:
    print("scorpion exists, skipping")
WASP = DB.find({"has_wasp": {"$gt": 0}})
try: 
    os.mkdir('wasp')
except:
    print("wasp exists, skipping")
# generate files for each species
print("Generating files for SPIDER...", end="")
gen_files(GENE_MAP, SPIDER, "spider")
print(Fore.GREEN + "DONE!" + Style.RESET_ALL)
print("Generating files for SCORPION...", end="")
gen_files(GENE_MAP, SCORPION, "scorpion")
print(Fore.GREEN + "DONE!" + Style.RESET_ALL)
print("Generating files for WASP...", end="")
gen_files(GENE_MAP, WASP, "wasp")
print(Fore.GREEN + "DONE!" + Style.RESET_ALL)


    

