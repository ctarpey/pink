### February 2016 
### Carolyn Tarpey 

####Python code to open a set of files in a folder that have the same extention
##and edit them by replacing the spaces with new lines, then writing them back out
### to a new file with the same name as the original, with an "edit" at the end. 

### Takes no arguments



import glob
path = "adjPvalues/*.txt"
for fname in glob.glob(path):
	nf = open(fname + "_edit.txt", 'w') 
	clean = open(fname).read().replace(' ', '\n')
	nf.write(str(clean))

