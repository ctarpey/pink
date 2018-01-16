#script convert a bayescan file to a bayenv2 file, reorder the allele counts
## July 29 2015
# carolyn tarpey 

## this needs two arguments:
#the first is the name of the bayescan file that has the allele counts, one loci per line, pops separated by [pop]=x (and the space replaced by tab)
#the second is the name of the output file you want (the bayenv2 file)


#!/bin/bash
import sys
import re

pat1 = r'[pop]='
AC_for_pop = dict()

with open(sys.argv[2], 'w') as OUTFILE:
	with open(sys.argv[1], 'r') as INFILE:
		#for i_ in range(4):
			#next(INFILE)
		for line in INFILE:
			if pat1 in line:
				popsplit = line.strip().split('=')
				current_pop = popsplit[1]
			elif line == '\n':
				continue
			else: 
				linesplit = line.strip().split("\t")
				current_locus = linesplit[1]
				allele_count_A = linesplit[4]
				allele_count_B = linesplit[5]
				#print current_pop,'\t', current_locus, '\t', allele_count_A
				#print current_pop,'\t', current_locus, '\t', allele_count_B
				#put the allele counts as the value in the dict with the locus as the key
				AC_of_locus = dict()
				AC_of_locus[current_locus] = (allele_count_A , allele_count_B)
				#put that dict in the overall dict with the pop as the key
				AC_for_pop[current_pop] = current_locus
				#OUTFILE.write('{}\t{}\t{}'.format(current_pop, current_locus, allele_count_A))
				#OUTFILE.write('{}\t{}\t{}'.format(current_pop, current_locus, allele_count_B))
			
			
			
			#for locus in loci:
				# allele_count_A = get_counts(pop, locus, allele)
				# allele_count_B = get_counts(pop, locus, allele)
				# AC_of_locus[locus] = (allele_count_A, allele_count_B)
			# AC_for_pop[pop] = AC_of_locus

count_of_A_pop1_locus3, count_of_B_pop1_locus3 = AC_for_pop[pop1][locus3] 


snplist_file.close()
outfile.close()