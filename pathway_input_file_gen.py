""" This script was created to search for EC reference in Uniprot,
    gather this info with others in Inovatoxin database and generate files.

    THe files generated are the input files for Pathway Tools.
"""
import re
import os
import datetime
import time
import pickle

# Retrieve Protein Names
from Bio import ExPASy, SeqIO
# Connect to MongoDB
from pymongo import MongoClient
# Nice Warnings
from colorama import Fore, Back, Style
# Progress bar
from tqdm import tqdm

class GeneNotFound(Exception):
    """ My own exception, for a gene that wasn't found."""
    pass

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
            print("\nGENE NOT FOUND. RETRYING (%d)" % retry )
            return get_protein_EC(gene, retry+1)
        return e
    except KeyboardInterrupt as k:
        print("\nKeyBoard Interrupt Signal received. Aborting")
        return k
    return None

def connect_Inovatoxin():
    """ Connection to Inovatoxin MongoDB """
    client = MongoClient()
    mongo = client['inovatoxin']
    return mongo['Proteins_protein']

def gen_files(gene_map, docs, species):
    """ Generates input files for PathwayTools, per species.
        Files are 3:
        1. genetic-elements.dat
            "Metadata" file, with all entries for a species.
            Points to the other files for each entry.
        2. fastas/ files
            Contains the sequence of that protein.
        3. infos/ files
            Contains addicional information of each entry.

        Only genes with an EC are added as an entry.
    """
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
    for doc in tqdm(docs):
        ec = gene_map.get(doc["Blast_fullAccession"], None)
        if ec is None:
            print(Fore.YELLOW + " EC NOT FOUND" + Style.RESET_ALL)
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

def collect_genes(genes, gene_map):
    """ Queries all genes missing and saves them in gene_map.
        The gene_map is: (gene_name, EC).
        If the gene has no EC, then it's (gene_name, None).

        Highly dependent on get_protein_EC function. If Exception is returned,
        then forwards the exception and terminates.
    """
    i = len(gene_map.keys())
    initial_time = datetime.datetime.now()
    print("CURRENT GENE -  TIME ELAPSED  -  AVERAGE TIME  - CURRENT GENE")
    for g in genes:
        time_before = datetime.datetime.now()
        i += 1
        print("\r [{:04}/{:04}] - {} - {} - Getting EC for {:12s} ".format(i, len(genes),
            time_before - initial_time, ((time_before - initial_time)/i), g), end="")
        ec = get_protein_EC(g)
        if isinstance(ec, Exception):
            return (gene_map, GeneNotFound(g + "not found."))
        else:
            gene_map[g] = ec
    print(Fore.GREEN + "All ECs collected!" + Style.RESET_ALL)
    return (gene_map, None)


def main():
    # Connect to DB
    DB = connect_Inovatoxin()
    # Get all names to search
    GENES = DB.distinct("Blast_fullAccession")
    # Gene cache 
    GENE_MAP = pickle.load(open("_gene_cache.pkl", "rb"))
    # removes already-searched genes from names to search
    for g in GENE_MAP.keys():
        GENES.remove(g)
    # search remaining names
    GENE_MAP, e = collect_genes(GENES, GENE_MAP)
    if e is not None:
        print("All collected proteins will be dumped in _gene_cache.pkl")
        pickle.dump(GENE_MAP, open("_gene_cache.pkl", "wb"))
        exit(1)
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

if __name__ == "__main__":
    main()

