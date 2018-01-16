#This is a script to take the genepop files that are made with iterative runs of populations and edit them to be used with Garrett's R code for filtering. 

import io
import re
import sys

f1 = open(sys.argv[1], "r") #open the genepop output file for reading
f1_lines = f1.readlines()
f1.close()

f2 = open(sys.argv[2], "w+") #this opens your output file and it creates the file if it doesn't exist

f2.write("Locus" + "\t" + "Fst" + "\n")

for line in f1_lines:
	if 'Locus:' in line:
		temp1= line.strip().split()
		#f2.write(temp1[2])
		print temp1
		#print temp1[2] #split line and write the second element
	if re.match(r'2      ', line): #if the line contains a sample, 
		temp2 = line.strip().split()
		#f2.write(temp2[2])
		#print temp2[2]
		#write it to the output file unchanged
		#f2.write("\t" + re.sub(r',', "\t", line))
	else: 
		continue
	#if it is the second line, (no pops or population name or stacks) write a tab, then split the remainder by replacing the , with tabs

f2.close()
