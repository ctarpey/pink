#!/usr/bin/perl -w
use strict;

my$MSTfile=$ARGV[0];
#my$CRIfile=$ARGV[1];

my%loci;
my$currentLG;
my$store_loci=0;
open(MST, "<$MSTfile")||die "cannot open $MSTfile:$!";
while(my $line=<MST>){
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
		$loci{$locus}=$currentLG;
		print "$locus\t$position\t$currentLG\n";
	}
}
close MST;
