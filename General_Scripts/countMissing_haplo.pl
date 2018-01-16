#!/usr/bin/perl -w
use strict;

my$header=<>;
chomp $header;
my($col1,$col2,$samples)=split "\t", $header, 3;
my%missing;
my@samples=split "\t", $samples;
while (my$line=<>){
	chomp $line;
	my($tag,$genoNum,$rest)=split "\t", $line, 3;
	my@genos=split "\t", $rest;
	foreach my $i (0..$#samples){
		if($genos[$i] eq "-"){
			$missing{$samples[$i]}++;
		}
	}
}

foreach my$key (keys %missing){
	print "$key\t$missing{$key}\n";
}