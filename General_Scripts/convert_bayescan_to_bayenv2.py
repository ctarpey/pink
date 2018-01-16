###Python code to convert a bayescan input file to Bayesenv2
###Carolyn Tarpey and Ryan Waples | September 2015

### typical workflow: genepop file, converted to GESTE/bayescan format with pgd spider, then manually edited
## delete the [loci] line at top and [populations] and find/replace spaces with tabs. 

###takes two arguments: 1) bayescan (GESTE format) input file, 2) output file name for Bayesenv2
# ie: convert_bayescan_to_bayenv2.py 16681_80_bayescan.txt 16681_80_bayenv2_in.txt

#!/bin/bash
import sys
import re

pat1 = r'[pop]='
pat2 = r'[populations]='
pat3 = r'[loci]='

AC_of_pop = dict()

with open(sys.argv[1], 'r') as INFILE:
	for line in INFILE:
		if pat1 in line:
			popsplit = line.strip().split('=')
			current_pop = int(popsplit[1])
			#print current_pop
			AC_of_locus = dict()
			
		elif line == '\n':
			continue
		elif pat2 in line:
			continue
		elif pat3 in line:
			continue
		else: 
			linesplit = line.strip().split("\t")
			current_locus = int(linesplit[0])
			allele_count_A = linesplit[3]
			allele_count_B = linesplit[4]
			#put the allele counts as the value in the dict with the locus as the key
			AC_list = [allele_count_A, allele_count_B]
			AC_of_locus[current_locus] = (AC_list)
			#put that dict in the overall dict with the pop as the key
			AC_of_pop[current_pop] = AC_of_locus

keys = AC_of_pop.keys()
print keys

with open(sys.argv[2], 'w') as OUTFILE:
	for locus in sorted(AC_of_locus.keys()):
		top_line = [] 
		bottom_line = [] 
		for pop in sorted(AC_of_pop.keys()):
			top_line.append(AC_of_pop[pop][locus][0])
			bottom_line.append(AC_of_pop[pop][locus][1])
		
		OUTFILE.write('\t'.join(top_line))
		OUTFILE.write('\t\n')
		OUTFILE.write('\t'.join(bottom_line))
		OUTFILE.write('\t\n')





############# test the order of the output data

# with open(sys.argv[2], 'w') as OUTFILE:
	# for locus in sorted(AC_of_locus.keys()):
		# top_line = [] 
		# bottom_line = [] 
		# for pop in sorted(AC_of_pop.keys()): 
			# top_line.append(str(pop))
			# top_line.append(str(locus))
			# bottom_line.append(str(pop))
			# bottom_line.append(str(locus))
		
		# OUTFILE.write('\t'.join(top_line))
		# OUTFILE.write('\t\n')
		# OUTFILE.write('\t'.join(bottom_line))
		# OUTFILE.write('\t\n')

# with open(sys.argv[2], 'w') as OUTFILE:
	# for locus in sorted(AC_of_locus.keys()):
		# top_line = [] 
		# bottom_line = [] 
		# for pop in sorted(AC_of_pop.keys()):
			# top_line.append(AC_of_pop[pop][locus][0])
			# top_line.append(str(locus))
			# bottom_line.append(AC_of_pop[pop][locus][1])
			# bottom_line.append(str(locus))
		
		# OUTFILE.write('\t'.join(top_line))
		# OUTFILE.write('\t\n')
		# OUTFILE.write('\t'.join(bottom_line))
		# OUTFILE.write('\t\n')					
		
#for pop, ac_dict in AC_of_pop.items():
#	print pop
#	print ac_dict.items()
#	print '\n'