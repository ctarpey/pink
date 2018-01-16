#!/usr/bin/perl -w
use strict;


my%loci;
my%LGlociNum;
my%LGlength;
my%LGbins;
my$lociNum=0;
my%binNum;
my$currentLG;
my$currentPosition;
my$store_loci=0;
while(my$line=<>){
	chomp $line;
	if($line=~/^group/){
		my($discard, $LG)=split " ", $line, 2;
		$currentLG=$LG;
	}
	if($line=~/\;BEGINOFGROUP/){
		$store_loci=1;
		next;
	}
	if($line=~/\;ENDOFGROUP/){
		$LGlociNum{$currentLG}=$lociNum;
		$LGlength{$currentLG}=$currentPosition;
		my$binNum=0;
		foreach my$key (keys %binNum){
			$binNum++;
		}
		$LGbins{$currentLG}=$binNum;
		%binNum=();
		$store_loci=0;
		$lociNum=0;
		$currentLG="";
		$currentPosition="";
	}
	if($store_loci==1){
		my($locus, $position)=split "\t", $line, 2;
		$loci{$locus}=$currentLG;
		$currentPosition=$position;
		$binNum{$position}++;
		$lociNum++;
	}
}

print "LinkageGroup\tLength(cM)\t#Markers\t#Bins\n";
foreach my$key (sort keys %LGlength){
	print "$key\t$LGlength{$key}\t$LGlociNum{$key}\t$LGbins{$key}\n";
}