### Plot the results from Lositan
###   Outlier tests of the population data
### Carolyn Tarpey | September 2015 
### ---------------------------------------



#
#Bare bones script to read LOSITAN data into R
#
# USAGE:
# The script expects to find:
#    a file called loci with locus He/Fst and
#    a file called ci with the confidence intervals
#
# A PNG graphic will be generated called output.png
#
# The graphic and the tables generated are 'bare bones' in
# the sense that they are provided as a STARTING POINT on
# how to read the data into R.
# User customization of the script is expected.
# Feel free to change and use it at your discretion
#
#(C) 2008 Tiago Antao
#This script is free software under the GPL v3
plot_ci <- function(f_name, bcolor, mcolor, tcolor) {
  cpl <- read.table(f_name, header=TRUE)
  lines(cpl[,1],cpl[,2], type='l', col=bcolor)
  lines(cpl[,1],cpl[,3], type='l', col=mcolor)
  lines(cpl[,1],cpl[,4], type='l', col=tcolor)
}

plot_loci <- function(f_name, color) {
  cpl <- read.table(f_name, header=TRUE, sep='\t')
  points(cpl[,2],cpl[,3], col=color)
}


png('output.png')
plot(-10,ylim=c(0,0.4),xlim=c(0,1), xlab='He', ylab='Fst')
plot_ci('ci', 'green', 'black', 'red')
plot_loci('loci', 'blue')
dev.off()
print('The purpose of this script is to provide a Bare bones example on how to read LOSITAN data')
print('You should customize it to your needs')
