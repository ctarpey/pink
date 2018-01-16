##Garrett wrote this, it takes the output of Genepop, an INF file formats it so you can gets the allele frequencies (I have never used it yet) 

#/usr/bin/perl -w
use strict;

my$keepSection=0;
my$rownum=0;
my$colnum=0;
my@rows;
my$locus;
my%popList;
print "Locus";
while(my$line=<>){
	chomp $line;
	if($line=~/Pop:/){
		my@columns=split " ", $line;
		my$pop=$columns[1];
		my($ID,$discard)=split "_", $pop;
		if(exists $popList{$ID}){
			next;
		}else{
			print "\t$ID";
			$popList{$ID}++;
		}
	}
	if(($line=~/^ Locus/)||($line=~/^Normal ending/)){
		if($keepSection==0){
			print "\n";
		}
		$line=~s/ //g;
		$colnum-=1;
		$rownum-=1;
		#print "$colnum\t$rownum\n";
		foreach my$i (0..$colnum){
			print "$locus";
			foreach my$j (0..$rownum){
				print "\t$rows[$j][$i]";
			}
			print "\n";
		}
		@rows=();
		my$discard;
		($discard,$locus)=split ":", $line;
		#print "$locus\n";
		$keepSection=1;
		$rownum=0;
		$colnum=0;
	}elsif($keepSection==1){
		if(($line=~/^ /)||($line=~/^-/)||($line=~/^$/)){
			next;
		}else{
			my@columns=split " ", $line;
			shift@columns;
			pop@columns;
			$colnum=scalar@columns;
			$rownum++;
			push@rows,\@columns;
			#foreach my$i (@columns){
			#	print "$i\t";
			#}
			#print "\n";
		}
	}
}