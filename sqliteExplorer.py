import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
from sqlite3 import Error

import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
 
 
def plotBarGraph(labels, results, run):
    """ Plots horizontal bar graph of labels vs results """
    y_pos = np.arange(len(labels))
    plt.barh(y_pos, results, align='center', alpha=0.5)
    plt.yticks(y_pos, labels)
    plt.xlabel('Quantidades Encontradas')
    plt.title("GO's de toxinas anotadas para {}".format(run))
    plt.show()


def clearData(results, labels):
    """ Removes empty values from lists, returning to clean lists"""
    # find empty values and mark them
    for i in range(len(results)):
        if results[i].size == 0:
            results[i] = 0
            labels[i] = ""
        else:
            results[i] = results[i][0]

    # remove marked values from both lists
    return(list(filter(lambda a: a != 0, results)), list(filter(lambda a: a != "", labels)))

def connect(filename):
    """ Connects to SQLite DB specified in filename."""
    try :
        conn = sqlite3.connect(filename)
    except Error as e:
        print(e)
        exit(1)
    return conn

def showTables(conn):
    """ Displays all tables in database connected by conn and some of the data in each."""
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    print("Available tables in DB")
    for t in tables:
        t = t[0]
        print("Tabela {}\n".format(t))
        df = pd.read_sql_query("select * from {} limit 5;".format(t), conn)
        print("Columns: \n")
        print(df.head())
        print("="*50 + "\n")

File = '/home/gguidini/Downloads/Trinotate.sqlite' # input("Path to sqlite file\n")
GOs = ['0046872', '0004888', '0019870', '0008200', '0008061', '0003950', '0005509', '0004568', '0030550', '0004872', '0019855', '0008092', '0005525', '0016881', '0004622', '0005515']
conn = connect(File)

# Get the table
df = pd.read_sql_query("SELECT DISTINCT * FROM final;", conn)
results = []
labels = []
# extract index of rows that contain targets
for go in GOs:
    results.append(df[df['id'] == ('GO:'+go)]['id'].value_counts().values)
    labels.append('GO:'+go)

(results, labels) = clearData(results, labels)
plotBarGraph(labels, results, "Escorpião, par 1")










