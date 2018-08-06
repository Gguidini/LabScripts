"""
To quickly filter possible toxins by Full Recname from Blast hits among various hits.
Basicly searches for some words in a .xls file. Results moved to out file.
"""

import re
import os

file = input("Input name\n")
toxin = re.compile(r'toxin', re.IGNORECASE) 
poison = re.compile(r'poison', re.IGNORECASE)
venom = re.compile(r'venom', re.IGNORECASE)
out = open(input("Output name\n"), 'w')

with open(file, 'r') as fd:
	line = fd.readline()
	out.write(line) # primeira linha sempre existe 
	cnt = 1
	while line:
		booly = False
		line = fd.readline()
		if toxin.search(line):
			booly = True
		elif poison.search(line):
			booly = True
		elif venom.search(line):
			booly = True

		if booly:
			out.write(line)
			cnt += 1
	print("Found {} matching lines".format(cnt))
