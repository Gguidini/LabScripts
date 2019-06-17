import re
import sqlite3
import sys

if __name__ == "__main__":
    # Connects to database
    if len(sys.argv) > 1:
        FILE = sys.argv[1]
    else:
        FILE = input("Input file (SignalP summary):")

    # DATABASE = input("Database file:")
    # CONN = sqlite3.connect(DATABASE)
    # CURSOR = CONN.cursor()
    # # Remove all entries in Database and Recreate the table
    # CREATE = """
    # DROP TABLE IF EXISTS SignalP;
    # CREATE TABLE SignalP(query_prot_id,start REAL,end REAL,score REAL,prediction, CS_type);
    # CREATE UNIQUE INDEX QueryID ON SignalP(query_prot_id);
    # """
    # CURSOR.executemany(CREATE)
    # CONN.commit()
    # Creates new entries
    with open(FILE, "r") as fd:
        for line in fd:
            if line.startswith("#"):
                continue
            parts = line.split('\t')
            if parts[1] == 'OTHER':
                query = """
                INSERT INTO
                    SignalP(query_prot_id, start, end, score, prediction)
                VALUES
                    (?, NULL, NULL, ?, 'OTHER');
                """
                args = (parts[0], float(parts[2]))
            else:
                query = """
                INSERT INTO
                    SignalP(query_prot_id, start, end, score, prediction)
                VALUES
                    (?, ?, ?, ?, 'SP(Sec/SPI)');
                """
                args = (parts[0], float(parts[2]))
