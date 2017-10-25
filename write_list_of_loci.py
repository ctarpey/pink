#script to write the number of loci that we have in pink catalog 4 as a list, one loci per line
# open the output file and make sure its in Linux EOL, and delete the first line. 

loci_list=list(range(997122))

print '\n'.join(str(p) for p in loci_list)

