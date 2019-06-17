import re
import sqlite3
import sys

if __name__ == "__main__":
    # Connects to database
    if len(sys.argv) > 1:
        FILE = sys.argv[1]
    else:
        FILE = input("Input file (SignalP summary):")

    DATABASE = input("Database file:")
    CONN = sqlite3.connect(DATABASE)
    CURSOR = CONN.cursor()
    # Remove all entries in Database and Recreate the table
    CREATE = [
    'DROP TABLE IF EXISTS SignalP;',
    """CREATE TABLE SignalP(
        query_prot_id TEXT,
        start INTEGER,
        end INTEGER,
        score REAL,
        prediction TEXT,
        CS_type TEXT,
        PRIMARY KEY(query_prot_id));""",
    'CREATE UNIQUE INDEX QueryID ON SignalP(query_prot_id);'
    ]
    for q in CREATE:
        CURSOR.execute(q)
    CONN.commit()
    RGX = re.compile(r"CS pos: (\d+)-(\d+)\. ([\w-]+)\.")
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
                args = (parts[0], float(parts[3]))
            else:
                query = """
                INSERT INTO
                    SignalP(query_prot_id, start, end, score, prediction, CS_type)
                VALUES
                    (?, ?, ?, ?, 'SP(Sec/SPI)', ?);
                """
                m = re.match(RGX, parts[4]).groups()
                args = (parts[0], int(m[0]), int(m[1]), float(parts[2]), m[2])
            CURSOR.execute(query, args)
        CONN.commit()
    CONN.close()
