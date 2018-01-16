#This is a script to take one genepop file that does not have one locus per line, and edit the formatting so that it can be imported to R as a table. 
#this runs with python, name of the script, genepop file to edit and the desired output file name

import io
import re
import sys


#f2 = sys.argv[2] Output file
#f1 = sys.argv[1] Genepop file to edit
f2 = open(sys.argv[2], "w+")
with io.open(sys.argv[1], 'r') as f1:
	for line in f1:
		if 'Stacks' in line:
			continue #if line contains('Stacks') skip it
		if re.match(r'pop', line): #if the line contains pop, skip it
			continue
		if re.match(r'P', line): #if the line contains a sample, 
			f2.write(line) #write it to the output file unchanged
		else: 
			f2.write("\t" + re.sub(r',', "\t", line))
		#if it is the second line, (no pops or population name or stacks) write a tab, then split the remainder by replacing the , with tabs
	f1.close()
	f2.close()
