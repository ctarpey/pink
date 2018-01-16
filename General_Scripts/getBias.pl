#!usr/bin/perl -w

use strict;

my$locus=42841;
my$outputFile="alleleBias.txt";
my@files=`ls`;

my%whitelist;
while(my$line=<>){
	chomp $line;
	$whitelist{$line}++;
	#print "$line\n";
}
#store genotypes for loci in whitelist
my$haplotypeFile="batch_10.haplotypes.tsv";
open(HAPLOTYPES,"<$haplotypeFile")||die "cannot open $haplotypeFile:$!";
my$header=<HAPLOTYPES>;
chomp $header;
my@samples=split "\t", $header;
my%genotypes;
while(my$line=<HAPLOTYPES>){
	chomp $line;
	my%individualGenotypes;
	my@genotypes=split "\t", $line;
	my$marker=$genotypes[0];
	#print "$marker\n";
	if(exists $whitelist{$marker}){
		foreach my$i (0..$#genotypes){
			$individualGenotypes{$samples[$i]}=$genotypes[$i];
		}
		#print "$marker\n";
		$genotypes{$marker}=\%individualGenotypes;
	}
}
close HAPLOTYPES;
print "Genotypes read in\n";

foreach my$marker (sort keys %genotypes){
	#print "$marker\t";
	foreach my$individual (sort keys %{$genotypes{$marker}}){
		#print "$individual:\t$genotypes{$marker}{$individual}\t";
	}
	#print "\n";
}

#store alleles for each locus in whitelist
my%alleles;
foreach my$i (@files){
	my%locusAlleles;
	chomp $i;
	if($i=~/matches/){
		open(FILE,"<$i")||die "cannot open file $i:$!";
		while(my$line=<FILE>){
			chomp $line;
			my@columns=split "\t", $line,8;
			if(exists $whitelist{$columns[2]}){
				$locusAlleles{$columns[5]}++;
				#$alleles{$columns[2]}=\%locusAlleles;
				$alleles{$columns[2]}{$columns[5]}++;
			}
		}
		close FILE;
	}
}
print "Alleles discovered\n";

my%alleleCounts;
foreach my$marker (sort keys %alleles){
	#print "$marker\t";
	foreach my$variant (sort keys %{$alleles{$marker}}){
		#print "$variant\t";
		$alleleCounts{$marker}++;
	}
	#print "\n";
}

#open files and get biases for each locus in whitelist
my%markerInfo;
foreach my$i (@files){
	chomp $i;
	if($i=~/matches/){
		my$file=$i;
		my($individual,$rest)=split '\.', $file, 2;
		open(FILE,"<$i")||die "cannot open file $i:$!";
		while(my$line=<FILE>){
			chomp $line;
			my@columns=split "\t", $line,8;
			if(exists $whitelist{$columns[2]}){
				#allele is $columns[6], count is $columns[6]
				$markerInfo{$columns[2]}{$individual}{$columns[5]}=$columns[6];
			}
		}
	}
}
print "Allele Counts obtained\n";
#print markers, alleles, and counts for each individual
my%ratios;
foreach my$marker (sort keys %whitelist){
	#print "$marker\t";
	foreach my$individual (sort keys %{$markerInfo{$marker}}){
		#print "$individual\t";
		my@counts;
		my@variants;
		foreach my$variant (sort keys %{$alleles{$marker}}){
			if(exists $markerInfo{$marker}{$individual}{$variant}){
				#print "$variant\t$markerInfo{$marker}{$individual}{$variant}\t";
				push@counts,$markerInfo{$marker}{$individual}{$variant};
				push@variants,$variant;
			}else{
				$markerInfo{$marker}{$individual}{$variant}=0;
				#print "$variant\t$markerInfo{$marker}{$individual}{$variant}\t";
				push@counts,$markerInfo{$marker}{$individual}{$variant};
				push@variants,$variant;
			}
		}
		#my$variantNum=scalar@variants;
		#print "variants:$variantNum\t";
		foreach my$i (0..$#counts){
			#print "$i\t";
			foreach my$j (1..$#counts){
				if(($i==$j)||($i>$j)){
					next;
				}else{
					if(($counts[$i]>0)||($counts[$j]>0)){
						my$ratio;
						if($counts[$j]==0){
							$ratio=1;
						}else{
							$ratio=$counts[$i]/$counts[$j];
							my$comparison=$variants[$i]."/".$variants[$j];
							#print "$comparison\t";
							#print "$counts[$i]/$counts[$j]\t";
							if($genotypes{$marker}{$individual}=~/\//){
								if($counts[$i]>$counts[$j]){
									push@{$ratios{$marker}{$comparison}},1;
								}else{
									push@{$ratios{$marker}{$comparison}},0;
								}
							}
						}
					}
				}
			}
		}
		#print "\n";
	}
	#print "\n";
}

foreach my$marker (sort keys %whitelist){
	#print "$marker\t";
	foreach my$comparison (sort keys %{$ratios{$marker}}){
		#print "$comparison\t";
		my$sum=0;
		my$count=0;
		foreach my$i (0..$#{$ratios{$marker}{$comparison}}){
			$sum+=$ratios{$marker}{$comparison}[$i];
			$count++;
			#print "$ratios{$marker}{$comparison}[$i]\t";
		}
		my$ratio=$sum/$count;
		my$rounded = sprintf "%.2f", $ratio;		
		#print "$sum/$count\t$rounded\t";
		print "$marker\t$comparison\t$sum/$count\t$rounded\n";
	}
	#print "\n";
}
