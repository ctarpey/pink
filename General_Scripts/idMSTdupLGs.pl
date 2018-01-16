#!/usr/bin/perl -w
use strict;

my%loci;
my$store_loci=0;
my$currentLG;
while(my$line=<>){
	chomp $line;
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
		$store_loci=0;
		$currentLG="";
	}
	if($store_loci==1){
		my($locus, $position)=split "\t", $line, 2;
		$locus=~s/R//;
		if(exists $loci{$locus}){
			$loci{$locus}.="\t".$currentLG;
		}else{
			$loci{$locus}=$currentLG;
		}
	}
}

my%groupings;
foreach my$key (keys %loci){
	$groupings{$loci{$key}}++;
}

foreach my$key (keys %groupings){
	print "$key\n";
}