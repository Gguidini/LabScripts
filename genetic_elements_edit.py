import itertools

PATH = 'MY_PATH/TO_FILES/'

def main():
    out_file = open('genetic-elements.ok.dat', 'w')
    with open('genetic-elements.dat', 'r') as infile:
        line = infile.readline()
        while line:
            if line.startswith('SEQ-FILE'):
                line = line[12:]
                line = 'SEQ-FILE    ' + PATH + line
            elif line.startswith('ANNOT-FILE'):
                line = line[14:]
                line = 'ANNOT-FILE    ' + PATH + line
            out_file.write(line)
            line = infile.readline()
    out_file.close()

if __name__ == "__main__":
    main()
