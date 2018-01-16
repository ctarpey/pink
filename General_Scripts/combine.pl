#/usr/bin/perl -w


use strict;

my@files=`ls PSBC1`;
my%samples;
foreach my$file (@files){
	if($file=~/.fq.gz/){
		my@info=split '\.', $file;
		#print "$info[0]\n";
		$samples{$info[0]}++;
	}
}

foreach my$sample (sort keys %samples){
	print "$sample\n";
	my$catCommand1="cat ./PSBC1/".$sample.".1.fq ./PSBC1r/".$sample.".1.fq ./PSBC2/".$sample.".1.fq ./PSCB2r/".$sample.".1.fq >./Combined/".$sample.".1.combined.fq";
	print "$catCommand1\n";
	#`$catCommand1`;
	my$catCommand2="cat ./PSBC1/".$sample.".2.fq ./PSBC1r/".$sample.".2.fq ./PSBC2/".$sample.".2.fq ./PSBC2r/".$sample.".2.fq >./Combined/".$sample.".2.combined.fq";
	print "$catCommand2\n";
	#`$catCommand2`;
}