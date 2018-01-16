#!/usr/perl -w
use strict;

my$samples=0;
print "Sample\tGenotyped Loci\tMissing Loci\tPercent Genotyped\n";
while(my$line=<>){
	chomp $line;
	if($line eq "Pop"){
		$samples=1;
		next;
	}elsif($samples==1){
		my@columns=split "\t", $line;
		my$sample=$columns[0];
		my$missing=0;
		my$columns=scalar@columns;
		my$lociNum=$columns-1;
		foreach my$i (1..$#columns){
			if($columns[$i] eq "0000"){
				$missing++;
			}
		}
		my$genotyped=$lociNum-$missing;
		my$percGeno=100-($missing/$lociNum*100);
		print "$sample:\t$genotyped\t$missing\t$percGeno\n";
	}
}