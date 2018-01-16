#This is a script to take the genepop files that are made with iterative runs of populations and edit them to be used with Garrett's R code for filtering. 

import io
import re

for i in range(1,30): #change this to the # of whitelists you have
	#this opens your output file and it creates the file if it doesnt exist
	f2 = open("whitelist_" + str(i) + ".genepop", "w+")
	with io.open(r'whitelist_' + str(i) + r'\batch_4.genepop', 'r') as f1:
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
