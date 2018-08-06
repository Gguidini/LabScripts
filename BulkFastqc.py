"""
This program is meant to analyse all the .fastq files in a directory using FASTQC
The folder to be processed can be passed in the command call as argument
Otherwhise, it'll be asked to the user
RUN USING python3
by Gguidini, 2018
"""

import os
import subprocess as sub
import sys
import fnmatch
import time
from threading import Thread
import queue
import time

processes_running = 0
done = 0

class Th(Thread):
	def __init__ (self, fileName):
		Thread.__init__(self)
		self.f = fileName

	def run(self):
		global processes_running
		processes_running += 1
		op = sub.Popen(["fastqc", f], stdout=sub.PIPE)
		op.communicate()	# forces thread to wait until cmd is completed
		processes_running -= 1
		info(f)

#TODO: Information not showing correctly: name of finished programs gets repeated.
def info(fileName):
	""" Informs user of the progress of operations"""
	global done
	global l
	global filesDone

	if fileName != "":
		done += 1
		filesDone.append(fileName)
	sub.run("clear")    # clear screen
	print("[{0:^3}/{1:^3}] Files processed".format(len(filesDone), len(l)))
	print("Files finished:")
	for f in filesDone:
		print(f)

def purge(Names, Pattern):
    """ Removes all files from current folder that match pattern"""
    for f in Names:
        if fnmatch.fnmatch(f, Pattern):
            os.remove(f)

# Path to directory
path = ""
if len(sys.argv) == 1:
    path = input("What folder do you want to FASTC?")
else:
    path = sys.argv[1]

# Goes to specified folder
go = False
while not go:
    try:
        os.chdir(path)
        go = True
    except (FileNotFoundError, NotADirectoryError):
        print("The directory was not found. Directories in this folder are:")
        with os.scandir(path) as it:
            for entry in it:
                if not entry.name.startswith('.') and entry.is_dir():
                    print(entry.name)
        path = input("Choose another directory")

# Proceed of not
print("Ok, I'm the correct directory now:")
print(os.getcwd())
r = input("Proceed? (Y/n)")
if r == 'n' or r == 'N':
    print("Aborting... bye")
    exit()

# Get list of files in directory
l = os.listdir()
# filters list so only .fastq will be there
l = fnmatch.filter(l, "*.fastq")
done = 0
# queue of processes
Q = queue.Queue(len(l)+1)
for f in l:
	Q.put_nowait(f) # filling the queue
# runs commands
filesDone = []
while not Q.empty():
	if processes_running < 3:
		f = Q.get_nowait()
		t = Th(f)    
		t.start()
		time.sleep(2) # famous gambiarra
# removes the .zip created on the process
# purge(os.listdir(), "*.zip")





