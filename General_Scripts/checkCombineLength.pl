#/usr/bin/perl -w
use strict;

my@files=`ls PSBC1`;
my%samples;
foreach my$file (@files){
	if($file=~/.fq/){
		my@info=split '\.', $file;
		$samples{$info[0]}++;
	}
}

print "Sample\tPSBC1\tPSBC1r\tPSBC2\tPSBC2r\tCombined\tMatch\n";
foreach my$sample (sort keys %samples){
	my$lengthCommand1="wc -l ./PSBC1/".$sample.".1.fq";
	my$lengthCommand2="wc -l ./PSBC1r/".$sample.".1.fq";
	my$lengthCommand3="wc -l ./PSBC2/".$sample.".1.fq";
	my$lengthCommand4="wc -l ./PSBC2r/".$sample.".1.fq";
	my$lengthCommand5="wc -l ./combined/".$sample.".1.combined.fq";
	my$lengthResult1=`$lengthCommand1`;
	my$lengthResult2=`$lengthCommand2`;
	my$lengthResult3=`$lengthCommand3`;
	my$lengthResult4=`$lengthCommand4`;
	my$lengthResult5=`$lengthCommand5`;
	chomp $lengthResult1;
	chomp $lengthResult2;
	chomp $lengthResult3;
	chomp $lengthResult4;
	chomp $lengthResult5;
	my$match1;
	if($lengthResult5==$lengthResult1+$lengthResult2+$lengthResult3+$lengthResult4){
		$match1=1;
	}else{
		$match1=0;
	}
	print "$sample.1.fq\t$lengthResult1\t$lengthResult2\t$lengthResult3\t$lengthResult4\t$lengthResult5\t$match1\n";
	
	my$lengthCommand6="wc -l ./PSBC1/".$sample.".2.fq";
	my$lengthCommand7="wc -l ./PSBC1r/".$sample.".2.fq";
	my$lengthCommand8="wc -l ./PSBC2/".$sample.".2.fq";
	my$lengthCommand9="wc -l ./PSBC2r/".$sample.".2.fq";
	my$lengthCommand10="wc -l ./combined/".$sample.".2.combined.fq";
	my$lengthResult6=`$lengthCommand6`;
	my$lengthResult7=`$lengthCommand7`;
	my$lengthResult8=`$lengthCommand8`;
	my$lengthResult9=`$lengthCommand9`;
	my$lengthResult10=`$lengthCommand10`;
	chomp $lengthResult6;
	chomp $lengthResult7;
	chomp $lengthResult8;
	chomp $lengthResult9;
	chomp $lengthResult10;
	my$match2;
	if($lengthResult10==$lengthResult6+$lengthResult7+$lengthResult8+$lengthResult9){
		$match2=1;
	}else{
		$match2=0;
	}
	print "$sample.2.fq\t$lengthResult6\t$lengthResult7\t$lengthResult8\t$lengthResult9\t$lengthResult10\t$match2\n";
}