###Python code to edit the PGD ped file to the LFMM ped file
###Carolyn Tarpey and Ryan Waples | September 2015

#workflow: get a genepop file and convert it to ped with pgd spider and use that as the input file here
###takes two arguments: 1) ped input file from PGDspider, 2) output file name for LFMM ped file
# ie: edit_PGDped_to_LFMMped.py LFMM_ped_pgd.ped LFMM_ped.txt

## Here is an example of what the output format should be: 
#1 	SAMPLE0 0 0 2 2 1 2 3 3 1 1 2 1
#2 	SAMPLE1 0 0 1 2 2 1 1 3 0 4 1 1
#3 	SAMPLE2 0 0 2 1 2 2 3 3 1 4 1 2

#!/bin/bash
import sys
import re

pat1 = r'pop_1 '
pat2 = r'pop_2 '
pat3 = r'pop_3 '
pat4 = r'pop_4 '
pat5 = r'pop_5 '
pat6 = r'pop_6 '
pat7 = r'pop_7 '
pat8 = r'pop_8 '
pat9 = r'pop_9 '
pat10 = r'pop_10 '
pat11 = r'pop_11 '
pat12 = r'pop_12 '
pat13 = r'pop_13 '
pat14 = r'pop_14 '

pat1a = r'AMUR10 \t '
pat2a = r'AMUR11 \t '
pat3a = r'HAYLY09 \t '
pat4a = r'HAYLY10 \t '
pat5a = r'KOPE91 \t '
pat6a = r'KOPE96 \t '
pat7a = r'KUSHI06 \t '
pat8a = r'KUSHI07 \t '
pat9a = r'NOME91 \t '
pat10a = r'NOME94 \t '
pat11a = r'SNOH03 \t '
pat12a = r'SNOH96 \t '
pat13a = r'TAUY09 \t '
pat14a = r'TAUY12 \t '

pat01 = r'01'
pat02 = r'02'
pat03 = r'03'
pat04 = r'04'

pat01a = r'1'
pat02a = r'2'
pat03a = r'3'
pat04a = r'4'


with open(sys.argv[1], 'r') as INFILE:
	with open(sys.argv[2], 'w') as OUTFILE:
		for line in INFILE:
			sample = line.strip()
			if pat1 in line:
				sample = re.sub(pat1, pat1a, sample)
			if pat2 in line:
				sample = re.sub(pat2, pat2a, sample) 
			if pat3 in line:
				sample = re.sub(pat3, pat3a, sample)
			if pat4 in line:
				sample = re.sub(pat4, pat4a, sample)
			if pat5 in line:
				sample = re.sub(pat5, pat5a, sample)
			if pat6 in line:
				sample = re.sub(pat6, pat6a, sample)
			if pat7 in line:
				sample = re.sub(pat7, pat7a, sample)
			if pat8 in line:
				sample = re.sub(pat8, pat8a, sample)
			if pat9 in line:
				sample = re.sub(pat9, pat9a, sample)
			if pat10 in line:
				sample = re.sub(pat10, pat10a, sample)
			if pat11 in line:
				sample = re.sub(pat11, pat11a, sample)
			if pat12 in line:
				sample = re.sub(pat12, pat12a, sample)
			if pat13 in line:
				sample = re.sub(pat13, pat13a, sample)
			if pat14 in line:
				sample = re.sub(pat14, pat14a, sample)
			sample = sample.split(" ")
			frontslice = sample[0:6]
			cutfrontslice = frontslice[0:4]
			backslice = sample[6:]
			backslicedit1 = re.sub(pat01, pat01a, " ".join(backslice))
			backslicedit2 = re.sub(pat02, pat02a, backslicedit1)
			backslicedit3 = re.sub(pat03, pat03a, backslicedit2)
			backslicedit4 = re.sub(pat04, pat04a, backslicedit3)
			backslicedit4 = backslicedit4.split(" ")
			editedline = cutfrontslice + backslicedit4
			OUTFILE.write(" ".join(editedline))
			OUTFILE.write("\n")
			
