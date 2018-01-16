#!/usr/bin/perl -w
use strict;

my$header=<>;
chomp $header;
my@headers=split "\t", $header;
print "$header\n";

while(my$line=<>){
	my$consensus=0;
	my$other=0;
    chomp $line;
    my@genos=split "\t", $line;
    foreach my$i (2..$#genos){
        if ($genos[$i] =~ /consensus/){
			$consensus=1;
        }
    }
	unless($consensus==1){
		print "$line\n";
	}
}

