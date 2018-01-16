#/usr/bin/perl -w
use strict;

my$whitelistCount=1;
my$interval=10000;
my$outfile="whitelist_".$whitelistCount.".txt";
open(OUTFILE,">$outfile")||die "cannot open $outfile:$!";
my$count=0;
while(my$line=<>){
	chomp $line;
	$count++;
	print OUTFILE "$line\n";
	if($count%$interval==0){
		#print OUTFILE "$line\n";
		close OUTFILE;
		$whitelistCount++;
		$outfile="whitelist_".$whitelistCount.".txt";
		open(OUTFILE,">$outfile")||die "cannot open $outfile:$!";
		#$whitelistCount++;
		#print "$line\t$whitelistCount\n";
	}
	
}
close OUTFILE;