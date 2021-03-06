####python plan for converting table s2 to FASTA format
## feb 18 2015 
# carolyn tarpey (&garrett)

## this needs two arguments:
#the first is the name of the excel file that needs converting 
#the second is the name of the output file you want (the FASTA file)

#!/bin/bash

import sys
import re

#open file
excel_file = open(sys.argv[1], "r")
FASTA = open(sys.argv[2],"w") 

for line in excel_file:#read one line  of the excel file at a time and 
	columns = line.split("\t")#take that line and split it up by the tabs
	#print columns 
	newline =[ ">", "\t", columns[1], "\t", columns[2], "\t", columns[3], "\n" ] #> the second column tab third column tab fourth column end line
	#print newline
	FASTA.write(''.join(newline)) # write this to the output file: > the second column tab third column tab fourth column end line
	seq = columns[4]
	pat1 = r'(\[)'
	pat2 = r'(/\w])'
	seq_new = ""
	x = re.sub(pat1, seq_new, seq)
	#print x
	y = re.sub(pat2, seq_new, x)
	#print y
	
	FASTA.write(''.join(y)) # write this to the output file: the edited sequence and the end of line 
	FASTA.write("\n") #skip a line in the output

excel_file.close()
FASTA.close()


