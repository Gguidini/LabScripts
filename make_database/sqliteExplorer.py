import pandas as pd
import sqlite3
from sqlite3 import Error

import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#import seaborn as sns
def plotGraph(x, n):
    ax = x.plot(kind='barh', title='GOs de toxinas anotadas para {}'.format(n), figsize=(6, 3.5))
    plt.tight_layout()
    plt.legend(loc='best')
    return ax

def anotateNumbers(ax):
    # set individual bar lables using above list
    for i in ax.patches:
        # get_width pulls left or right; get_y pushes up or down
        ax.text(i.get_width()+ 1, i.get_y()+.1, \
                str(int((i.get_width()))), fontsize=8, color='dimgrey')

    # invert for largest on top 
    ax.invert_yaxis()

    # invert for largest on top 
    ax.invert_yaxis()

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

def getTables(conns, names, GOs, r):
    if names != None and len(conns) != len(names):
        print("Names and connections not the same length!")
        return None

    for i, c in enumerate(conns):
        df = pd.read_sql_query("SELECT id FROM final;", c)
        df = df['id'].value_counts().filter(items=['GO:'+go for go in GOs])
        print(df)
        print('='*150)
        r[names[i]] = df

    return r

def mergeColumn(df, index, new_name):
    x = df[index[0]]
    for i in index[1:]:
        x += df[i]
    df.drop(columns=index)
    df[new_name] = (x/3)
    return df


showTables(connect('/home/gguidini/Downloads/SQLite/Trinotate_scorpio_2.sqlite'))
exit(0)
# Files = input("Path to sqlite files (comma separated)\n").split(',')
# Files = list(map(str.strip, Files))
Files = ['/home/gguidini/Downloads/SQLite/Trinotate_scorpio_1.sqlite',
         '/home/gguidini/Downloads/SQLite/Trinotate_scorpio_2.sqlite',
         '/home/gguidini/Downloads/SQLite/Trinotate_scorpio_3.sqlite',
         '/home/gguidini/Downloads/SQLite/Trinotate_wasp_1.sqlite',
         '/home/gguidini/Downloads/SQLite/Trinotate_wasp_2.sqlite',
         '/home/gguidini/Downloads/SQLite/Trinotate_wasp_3.sqlite',
         '/home/gguidini/Downloads/SQLite/Trinotate_spider.sqlite']
# n = input(u"Qual a espécie?\n")
names = ['par 1', 'par 2', 'par 3']
GOs = ['0046872', '0004888', '0019870', '0008200', '0008061', '0003950', '0005509', '0004568', '0030550', '0004872', '0019855', '0008092', '0005525', '0016881', '0004622', '0005515']
conns = []
for f in Files:
    conns.append(connect(f))

#df = getTables(conns[0:3], names, GOs, pd.DataFrame(index=['GO:'+go for go in GOs]))
#df = mergeColumn(df, names, 'Escorpião')
#df = getTables(conns[3:6], names, GOs, df)
#df = mergeColumn(df, names, 'Vespa')
#df = getTables([conns[6]], ['Aranha'], GOs, df)
#df = df[['Escorpião', 'Aranha', 'Vespa']]
# Get the table

#sns.heatmap(df, annot=True, fmt='g')
#plt.title("Quantidades de GOs anotados por espécie")
#plt.tight_layout()
#plt.show()










