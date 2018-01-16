#This is a script combine VCF files. 
#this runs with python, name of the script, the file to add the lines to, the file with the lines you want to add to the other

import io
import re
import sys

#f1 = sys.argv[1] File to be appended to
#f2 = sys.argv[2] File to add to other

f1 = open(sys.argv[1], 'a')
with io.open(sys.argv[2], 'r') as f2:
	for line in f2:
		if '#CHROM' in line:
			continue #if line contains('#CHROM') skip it
		else: 
			f1.write(line)
f1.close()
f2.close()
