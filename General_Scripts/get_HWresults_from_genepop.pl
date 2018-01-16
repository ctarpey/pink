#/usr/bin/perl -w
use strict;

my$firstLocus=0;
my$getPval=0;
while(my$line=<>){
	chomp $line;
	if($line=~/^-/){
		next;
	}
	if($line=~/Locus/){
		$firstLocus++;
		my($discard,$locus)=split " ", $line, 2;
		$locus=~s/"//g;
		print "\n$locus";
	}
	if($getPval==1){
		my($pop,$pval,$rest)=split " ", $line, 3;
		print "\t$pval";
		#print "$pval";
	}
	if($line=~/POP/){
		$getPval=1;
	}
	if($line=~/^ *$/){
		$getPval=0;
		#print "\n";
	}
	#print "|";
}