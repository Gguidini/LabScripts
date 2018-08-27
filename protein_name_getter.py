"""
Script to move protein information from SQLite - created in the annotation process
to MongoDB - used in the website with Django.

Also unifies SQLite databases, extracts hit counts and retrieves protein names using BioPython.

WARNING: The script takes a lot of time to run depending on Entry size. TIme trackers over major 
time consuming functions
"""

# Retrieve Protein Names
from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO
# Connect to MongoDB
from pymongo import MongoClient
# Connect to SQLite
import sqlite3
from sqlite3 import Error
# Easily manage SQL Tables in transition
import pandas as pd
# Time tracking
import time
def getNameAndKeys(gene):
	""" 
	Gets name and keywords for Protein given gene.
	By Waldeyr Mendes, PhD
	"""
	try:
		with ExPASy.get_sprot_raw(gene) as handle:
			seq_record = SeqIO.read(handle, "swiss")
		fullName = seq_record.description[14:].split(";")[0]
		keys = seq_record.annotations["keywords"]
		return (fullName, keys)
	except :
		print("WARNING: %s not found" % gene)
		return (None, None)

def connect(filename):
    """ Connects to SQLite DB specified in filename."""
    try :
        conn = sqlite3.connect(filename)
    except Error as e:
        print(e)
        exit(1)
    return conn

def connectInovatoxin():
	""" Connection to MongoDB """
	client = MongoClient()
	db = client['inovatoxin']
	return db['Proteins_protein']

def uploadDocs(db, df, data):
	""" 
	Uploads docs to MongoDB in correct Collection.
	Docs are generated according to Django Model.
	"""
	print("Inserindo Documentos...")
	total = df.index
	columns = df.columns
	fields = ['Uniprot_accession', 'Uniprot_attribute_type', 'Pfam_accession', 'Pfam_domainName', 'Pfam_domainDescription', 'NCBI_taxonomyAccession', 'NCBI_taxonomyValue', 'Blast_percentIdentity', 'Blast_Evalue', 'Blast_fullAccession', 'Blast_GINumber', 'HMMER_domain', 'HMMER_domainDescription', 'HMMER_fullSeqEvalue', "transcript_sequence", 'orf_peptide', 'EggNOG_indexTerm', 'EggNOG_descriptionValue', 'GO_id', 'GO_name', 'GO_namespace', 'GO_def', 'has_scorpion', 'has_wasp', 'has_spider']
	for index, row in df.iterrows():
		print("[{}/{}] - {}\r".format(index, total, row['id']), end="")
		doc = {}
		t = data[row['FullAccession']]
		doc['id'] = index
		doc['name'] = t[0]
		doc['keywords'] = t[1]
		for idx, col in enumerate(columns):
			doc[fields[idx]] = row[col]
		# insert doc in DB
		db.insert_one(doc)
	print("Done!")


def getTable(conn):
	""" 
	Gets the necessary information from SQLite Databases into a pandas DataFrame
	so it can be manipulated.
	"""
	# insecure, but will be used locall, and only once.
	query = """ 
	SELECT DISTINCT UniprotIndex.Accession, UniprotIndex.AttributeType, PFAMreference.pfam_accession, PFAMreference.pfam_domainname, PFAMreference.pfam_domaindescription, TaxonomyIndex.NCBITaxonomyAccession, TaxonomyIndex.TaxonomyValue, BlastDbase.PercentIdentity, BlastDbase.Evalue, BlastDbase.FullAccession, BlastDbase.GINumber, HMMERDbase.HMMERDomain, HMMERDbase.HMMERTDomainDescription, HMMERDbase.FullSeqScore, Transcript.sequence, ORF.peptide, eggNOGIndex.eggNOGIndexTerm, eggNOGIndex.eggNOGDescriptionValue, go.id, go.name, go.namespace, go.def, UniprotIndex.LinkID
	FROM Transcript 
	INNER JOIN ORF ON ORF.transcript_id=Transcript.transcript_id
	INNER JOIN BlastDbase ON (BlastDbase.TrinityID=ORF.orf_id OR BlastDbase.TrinityID=ORF.transcript_id)
	INNER JOIN UniprotIndex ON  BlastDbase.FullAccession=UniprotIndex.Accession
	INNER JOIN go ON go.id=UniprotIndex.LinkID
	INNER JOIN HMMERDbase ON (HMMERDbase.QueryProtID=ORF.orf_id OR HMMERDbase.QueryProtID=ORF.transcript_id)
	INNER JOIN PFAMreference ON HMMERDbase.pfam_id=PFAMreference.pfam_accession
	LEFT JOIN TaxonomyIndex ON UniprotIndex.LinkId=TaxonomyIndex.NCBITaxonomyAccession
	LEFT JOIN eggNOGIndex ON eggNOGIndex.eggNOGIndexTerm=UniprotIndex.LinkId
	"""
	df = pd.read_sql_query( query, conn)
	return df

def makeCount(df):
	""" Returns hitcount for different Proteins' Ids"""
	return df['LinkId'].value_counts()


def getNames(accessions):
	""" Handles the process of retrieving names and keywords usign function getNameAndKey."""
	accessions = accessions.unique()
	print("Getting names...")
	total = accessions.size
	done = 0
	data = {}
	for acc in accessions:
		print("[{}/{}] - {}...\r".format(done, total, acc), end="")
		(f, k) = getNameAndKeys(acc)
		if k != None:
			k = ", ".join(k)
		#print(f, k)
		data[acc] = (f, k)
		done += 1
	print("[{}/{}] - Done!".format(done, total))
	return data

def getCountsTotal(counts):
	""" Adds hitcounts of same species """
	return counts[0].add(counts[1], fill_value=0).add(counts[2], fill_value=0).astype(int)

def main():
	t_inicial = time.time()
	prefix = input("Path to files? ")
	files = [ prefix + name for name in ['scorpio_1.sqlite', 'scorpio_2.sqlite', 'scorpio_3.sqlite', 'wasp_1.sqlite', 'wasp_2.sqlite', 'wasp_3.sqlite', 'spider.sqlite']]
	df = pd.DataFrame()
	counts = []

	done = 0
	print("Recuperando dados...")
	for f in files:
		print("[{}/3] - {}\r".format(done, f), end="")
		conn = connect(f)
		x = getTable(conn)
		counts.append(makeCount(x))
		df = pd.concat([df,x]).drop_duplicates().reset_index(drop=True)
		done += 1
	print("Tempo para recuperar informação do SQLite e juntar tabelas: {}".format(time.time() - t_inicial))

	# Adiciona contagem
	t_add = time.time()
	print("Analisando contagem de Ids")
	scorpio = getCountsTotal(counts[0:3])
	df['has_scorpion'] = [scorpio[df.iloc[idx]['LinkId']] for idx in df.index]
	wasp = getCountsTotal(counts[3:6])
	df['has_wasp'] = [wasp[df.iloc[idx]['LinkId']] for idx in df.index]
	spider = counts[6]
	df['has_spider'] = [spider[df.iloc[idx]['LinkId']] for idx in df.index]
	df = df.drop(columns=['LinkId'])
	print("Tempo de calculo de contagems (Ids): {}".format(time.time() - t_add))

	t_names = time.time()
	data = getNames(df['FullAccession'])
	print("Tempo de recuperar nomes: {}".format(time.time() - t_names))

	db = connectInovatoxin()
	t_upload = time.time()
	uploadDocs(db, df, data)
	print("Tempo de fazer upload: {}".format(time.time() - t_upload))

	print("Tempo total: {}".format(time.time() - t_inicial))

	#wasp = counts[3] + counts[4] + counts[5]
	#df['has_wasp'] = wasp
	#df['has_spider'] = counts[6]

if __name__ == '__main__':
	main()
