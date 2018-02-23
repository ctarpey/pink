####python plan for converting table s2 to FASTA format
## feb 18 2015 edited 20180222
# carolyn tarpey (&garrett)

## this needs two arguments:
#the first is the name of the excel file that needs converting 
#the second is the name of the output file you want (the FASTA file)

##CHANGE THE COLUMN INDEXES< and 0 is the first column, 1 is the second etc. 

#!/bin/bash

import sys
import re

#open file
excel_file = open(sys.argv[1], "r")
FASTA = open(sys.argv[2],"w") 

for line in excel_file:#read one line  of the excel file at a time and 
	columns = line.split("\t")#take that line and split it up by the tabs
	newline =[ ">", columns[0], "\n" ] #> the locus name
	#print newline
	FASTA.write(''.join(newline)) # write this to the output file: > the first column tab second column tab end line
	seq = columns[1]
	#print seq
	FASTA.write(''.join(seq)) # write this to the output file: the sequence
	FASTA.write("\n") #skip a line in the output

excel_file.close()
FASTA.close()


