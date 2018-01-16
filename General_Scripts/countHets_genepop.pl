#!/usr/perl -w
use strict;

my$samples=0;
print "Sample\tGenotyped Loci\tHeterozygous Loci\tPercent Heterozygous\n";
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
		my$het=0;
		foreach my$i (1..$#columns){
			my$allele1=substr $columns[$i], 0, 2;
			my$allele2=substr $columns[$i], 2, 2;
			if($allele1 ne $allele2){
				$het++;
			}
			#print "$allele1\t$allele2\n";
		}
		my$genotyped=$lociNum-$missing;
		my$percGeno=100-($missing/$lociNum*100);
		my$percHet=($het/$lociNum*100);
		my$percHetGeno=($het/$genotyped*100);
		#print "$sample:\t$genotyped\t$missing\t$percGeno\n";
		print "$sample:\t$genotyped\t$het\t$percHetGeno\n";
	}
}