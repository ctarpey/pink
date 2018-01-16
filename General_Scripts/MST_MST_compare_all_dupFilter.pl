#!/usr/bin/perl -w
use strict;

my$MSTfile1=$ARGV[0];
my$MSTfile2=$ARGV[1];
my$DuplicateFile=$ARGV[2];

my%LociDuplication;
open(DUPLICATION, "<$DuplicateFile")||die "cannot open $DuplicateFile:$!";
my$DupHeader=<DUPLICATION>;
chomp $DupHeader;
while(my$line=<DUPLICATION>){
	chomp $line;
	my($loci,$duplication)=split "\t", $line, 2;
	$LociDuplication{$loci}=$duplication;
}

my%MST1loci;
my$MST1currentLG;
my$MST1store_loci=0;
my%MST1LGlociNum;
my%MST1LGlength;
my$MST1lociNum=0;
my$MST1currentPosition;
open(MST1, "<$MSTfile1")||die "cannot open $MSTfile1:$!";
while(my $line=<MST1>){
	chomp $line;
	if($line=~/^group/){
		my($discard, $LG)=split " ", $line, 2;
		$MST1currentLG=$LG;
	}
	if($line=~/\;BEGINOFGROUP/){
		$MST1store_loci=1;
		next;
	}
	if($line=~/\;ENDOFGROUP/){
		$MST1LGlociNum{$MST1currentLG}=$MST1lociNum;
		$MST1LGlength{$MST1currentLG}=$MST1currentPosition;
		$MST1store_loci=0;
		$MST1currentLG="";
		$MST1currentPosition="";
		$MST1lociNum=0;
	}
	if($MST1store_loci==1){
		my($locus, $position)=split "\t", $line, 2;
		#$locus=~s/R//;
		if($LociDuplication{$locus}==1){
			$MST1loci{$locus}=$MST1currentLG;
			$MST1currentPosition=$position;
			$MST1lociNum++;
		}
	}
}
close MST1;

my%MST2loci;
my$MST2currentLG;
my$MST2store_loci=0;
my%MST2LGlociNum;
my%MST2LGlength;
my$MST2lociNum=0;
my$MST2currentPosition;
open(MST2, "<$MSTfile2")||die "cannot open $MSTfile2:$!";
while(my $line=<MST2>){
	chomp $line;
	if($line=~/^group/){
		my($discard, $LG)=split " ", $line, 2;
		$MST2currentLG=$LG;
	}
	if($line=~/\;BEGINOFGROUP/){
		$MST2store_loci=1;
		next;
	}
	if($line=~/\;ENDOFGROUP/){
		$MST2LGlociNum{$MST2currentLG}=$MST2lociNum;
		$MST2LGlength{$MST2currentLG}=$MST2currentPosition;
		$MST2store_loci=0;
		$MST2currentLG="";
		$MST2currentPosition="";
		$MST2lociNum=0;
	}
	if($MST2store_loci==1){
		my($locus, $position)=split "\t", $line, 2;
		#$locus=~s/R//;
		if($LociDuplication{$locus}==1){
			$MST2loci{$locus}=$MST2currentLG;
			$MST2currentPosition=$position;
			$MST2lociNum++;
			#print "$locus\t$MST2currentLG\n";
		}
	}
}
close MST2;

my%LGpairs;
foreach my$key (keys%MST1loci){
	my$Rkey=$key."R";
	if(((exists $MST1loci{$key})&&(exists $MST2loci{$key}))||((exists $MST1loci{$Rkey})&&(exists $MST2loci{$key}))){
		#print "$key\n";
		#print "$MST1loci{$key}\t$MST2loci{$key}\n";
		my$pair="$MST1loci{$key} $MST1LGlociNum{$MST1loci{$key}} $MST1LGlength{$MST1loci{$key}}\t$MST2loci{$key} $MST2LGlociNum{$MST2loci{$key}} $MST2LGlength{$MST2loci{$key}}";
		$LGpairs{$pair}++;
	}
}

foreach my$key (sort keys %LGpairs){
	print "$key\n";
}
