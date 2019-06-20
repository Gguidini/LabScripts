import sqlite3
import sys
import os
from tqdm import tqdm

def gen_files(species):
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
    # query to get extra info on enzime
    extra_info = """
    SELECT
        UniprotIndex.LinkId,
        UniprotIndex.AttributeType as 'Type'
    FROM
        UniprotIndex
    WHERE
        UniprotIndex.Accession = ? AND
        (UniprotIndex.AttributeType = 'D' OR UniprotIndex.AttributeType = 'G');
    """
    for doc in tqdm(ENZIMES, desc="Enzimas", unit="enz"):
        # Extract needed information for enzime
        (transcript_id, seq) = (doc[0], doc[4])
        if transcript_id is None:
            (transcript_id, seq) = (doc[1], doc[5])
        (acc, ec) = (doc[2], doc[3])
        CURSOR.execute(extra_info, (acc,))
        info = CURSOR.fetchall()
        go = None
        for row in info:
            if row[1] == 'D':
                name = row[0][14:]
            else:
                go = row[0]
        prot_id = transcript_id + "_" + species
        fasta_file = open(species + "/fastas/" + prot_id + ".fsa", 'w')
        info_file = open(species + "/infos/" + prot_id + ".pf", 'w')

        # General file
        text = "ID  " + prot_id + "\n"
        text += "NAME  " + name + "\n"
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
        text += "NAME  " + name + "\n"
        text += "PRODUCT-TYPE  P\n"
        text += "STARTBASE  1\n"
        text += "ENDBASE  " + str(len(seq)) + "\n"
        text += "EC  " + ec + "\n"
        if go is not None:
            text += "DBLINK  " + go + '\n'
        text += "//\n"
        info_file.write(text)
        info_file.close()
    genetic_elements.close()

if __name__ == "__main__":
    # Connects to database
    if len(sys.argv) > 1:
        DATABASE = sys.argv[1]
    else:
        DATABASE = input("Database file:")

    CONN = sqlite3.connect(DATABASE)
    CURSOR = CONN.cursor()
    SPECIES = input("What species is this? ")
    QUERY = """
    SELECT
        Transcript.transcript_id as 'Transcript',
        ORF.orf_id as 'ORF',
        BlastDbase.FullAccession as 'Accession',
        UniprotIndex.LinkId as 'EC Number',
        Transcript.sequence as 'Sequence',
        ORF.peptide as 'ORF'
    FROM
        UniprotIndex
    INNER JOIN BlastDbase ON
        BlastDbase.UniprotSearchString = UniprotIndex.Accession AND UniprotIndex.AttributeType = 'F'
    LEFT JOIN Transcript ON
        BlastDbase.TrinityID = Transcript.transcript_id
    LEFT JOIN ORF ON
        BlastDbase.TrinityID = ORF.orf_id
    GROUP BY 
        UniprotIndex.LinkId;
    """
    CURSOR.execute(QUERY)
    ENZIMES = CURSOR.fetchall()
    try:
        os.mkdir(SPECIES)
    except FileExistsError:
        print("Species directory already exists, skipping")
    gen_files(SPECIES)
    print("Completed.")