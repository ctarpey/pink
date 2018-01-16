#!/usr/bin/perl -w
use strict;

my$MSTfile=$ARGV[0];
my$markerList=$ARGV[1];

my%printMarkers;
open(MARKERLIST, "<$markerList")||die "cannot open $markerList:$!";
while(my$line=<MARKERLIST>){
	chomp $line;
	$printMarkers{$line}++;
}
close MARKERLIST;

open(MSTFILE, "<$MSTfile")||die "cannot open $MSTfile:$!";
my$lineCount=0;
while(my$line=<MSTFILE>){
	chomp $line;
	$lineCount++;
	my($locus,$rest)=split "\t", $line, 2;
	if($lineCount<=14){
		print "$line\n";
	}elsif(exists $printMarkers{$locus}){
		print "$line\n";
	}
}