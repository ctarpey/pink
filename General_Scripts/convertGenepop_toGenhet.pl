#/usr/bin/perl -w
use strict;

my$header=<>;
my$genoSection=0;
my$IDprint=0;
my@loci;
print "sampleid";
while(my$line=<>){
	chomp $line;
	if($line=~/pop/i){
		$genoSection=1;
		next;
	}
	if($genoSection==0){
		push@loci,$line;
	}else{
		if($IDprint==0){
			foreach my $i (@loci){
				my$allele1=$i."_1";
				my$allele2=$i."_2";
				print "\t$allele1\t$allele2";
			}
			print "\n";
			$IDprint=1;
		}
		my@columns=split "\t", $line;
		my$ID=$columns[0];
		$ID=~s/,//;
		print "$ID";
		foreach my$i (1..$#columns){
			my$allele1=substr $columns[$i],0,2;
			my$allele2=substr $columns[$i],2,2;
			if($allele1 eq "00"){
				$allele1="NA";
			}
			if($allele2 eq "00"){
				$allele2="NA";
			}
			print "\t$allele1\t$allele2";
		}
		print "\n";
	}
}