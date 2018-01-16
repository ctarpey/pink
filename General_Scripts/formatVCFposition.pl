#/usr/bin/perl -w
# runs with the vcf file as the input and a right arrow to the output file. 
use strict;

my$VCF=$ARGV[0];
chomp $VCF;

open(VCF,"<$VCF")||die "cannot open $VCF:$!";
my$locusCount=0;
my$lastTag="first";
while(my$line=<VCF>){
	chomp $line;
	if($line=~/^#/){
		print "$line\n";
	}else{
		my($chrom,$pos,$tag,$rest)=split "\t", $line,4;
		if($lastTag eq "first"){
			$lastTag=0;
		}elsif($tag!=$lastTag){
			$locusCount++;
			$lastTag=$tag;
		}
		my$correctedCount=$locusCount-1;
		my$correctedSNP=$pos-2;
		my$correctedPos=$correctedSNP-($correctedCount*94);
		#print "$tag\t$locusCount\t$pos\t$correctedPos\n";
		print "$line\n";
		#print "$chrom\t$correctedPos\t$tag\t$rest\n";
		print "$correctedPos\t$tag\n";
	}
}
close VCF;