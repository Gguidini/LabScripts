from pymongo import MongoClient

def connectInovatoxin():
    """ Connection to MongoDB """
    client = MongoClient()
    db = client['inovatoxin']
    return db['Proteins_protein']