#script to get union of multiple SNP lists
#default version of script works with 3 lists, will need to slightly modified for more files
#5/27/15
#wes larson, wlarson1@uw.edu

#define number of snp files
num_snp_files=3

#define file for each snplist, add more files and name them snplist4_file etc
snplist1_file="samp1.txt"
snplist2_file="samp2.txt"
snplist3_file="samp3.txt"

#make list of files, add new files to this list
snp_file_list=[snplist1_file,snplist2_file,snplist3_file]

#open output file, name the output file here
outfile=open("test_snplist_union.txt","w")

#define dictionary, will be snp name and num occurances
snp_dict={}
#iterate through each file, add to the dictionary
for i in snp_file_list:
    #open each file, make each line (snp) an entry in an array
    snplist_file=open(i,"r")
    snplist_array=snplist_file.readlines()
    snplist_file.close()
    #iterate through each line in the file, add to snp dict
    for j in snplist_array:
        current_snp=j.rstrip()
        if current_snp not in snp_dict.keys():
            snp_dict[current_snp]=1
        elif current_snp in snp_dict.keys():
            num_occur=snp_dict[current_snp]
            snp_dict[current_snp]=num_occur+1

for i in sorted(snp_dict.keys()):
    num_occ=snp_dict[i]
    if num_occ == num_snp_files:
        outfile.write(i+"\n")

outfile.close()