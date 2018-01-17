#pink_map_convert_to_FASTA_REMdup.py
#this was written by Ryan, and it makes a FASTA file from an excel sheet that has both the locus name and a sequence for that tag 


import os
import os.path
import numpy as np
import pandas as pd
from IPython.core.pylabtools import figsize
import matplotlib.pyplot as plt

# Try to read in this file G:\Analysis\Mapping\Align_to_chinook\PinkLinkagemap_consensus.xlsx
pink_ex = pd.read_excel('Z:\WORK\TARPEY\Exp_Pink_pops\Analysis\HWE_Alignment\Loci_FAILED_HWE.xlsx')
pink_ex.head(2)

#drop one of each duplicate
pink_ex.drop_duplicates(subset = ['Locus', 'Sequence'], inplace = True)

#check you only have one duplicate
len(set(pink_ex['Locus'])) == len(pink_ex['Locus'])

#True
### write out the consensus sequence from each locus as a FASTA file

for index, locus_row in pink_ex.iterrows():
    print '>{}\n{}'.format(locus_row['Locus'], locus_row['Sequence'])


