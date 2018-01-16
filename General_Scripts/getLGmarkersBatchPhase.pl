#!/usr/bin/perl -w
use strict;

my$MSTfile=$ARGV[0];
my$LGfile=$ARGV[1];

my%LGs;
open(LGFILE, "<$LGfile")||die "cannot open $LGfile:$!";
while(my$line=<LGFILE>){
	chomp $line;
	$LGs{$line}++;
}
close LGFILE;

open(MSTFILE, "$MSTfile")||die "cannot open $MSTfile:$!";
my%loci;
my$store_loci=0;
my$currentLG;
while(my$line=<MSTFILE>){
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
		#$locus=~s/R//;
		if(exists $LGs{$currentLG}){
			print "$locus\n";
		}
	}
}
