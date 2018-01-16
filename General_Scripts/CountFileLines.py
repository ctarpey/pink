#This is a script to check the combined length of any lengthy file. 
#this runs with python, name of the script, and the name of the file to count the lines of 

import io
import re
import sys

i = 0
with open(sys.argv[1]) as f:
	for i, l in enumerate(f):
		pass
	i + 1

print i

