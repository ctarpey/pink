#/usr/bin/perl -w

#To run it follow this example: perl countHetsMissing_genepop_byPop.pl inputfile >outputfile
# this does the calculation on a genepop input file, so the raw genotypes. 


use strict;

my$header=<>;

my@loci;
my%missingGenos;
my%hetGenos;
my%totalGenos;
my$locusSection=1;
my$pop=0;
print "Pop\tID\tTotalGenos\tMissingGenos\tHeterozygotes\tMissingPerc\tHetPerc\n";
while(my$line=<>){
	chomp $line;
	if($line=~/pop/i){
		$locusSection=0;
		if($pop==0){
			$pop++;
			next;
		}else{
			foreach my$key (sort keys %totalGenos){
				my$missingPerc;
				my$hetPerc;
				if(not exists $missingGenos{$key}){
					$missingGenos{$key}=0;
					$missingPerc=0
				}else{
					$missingPerc=$missingGenos{$key}/$totalGenos{$key};
				}
				if(not exists $hetGenos{$key}){
					$hetGenos{$key}=0;
					$hetPerc=0;
				}else{
					$hetPerc=$hetGenos{$key}/$totalGenos{$key};
				}
				print "$pop\t$key\t$totalGenos{$key}\t$missingGenos{$key}\t$hetGenos{$key}\t$missingPerc\t$hetPerc\n";
			}
			%missingGenos=();
			%hetGenos=();
			%totalGenos=();
			$pop++;
			next;
		}
	}
	if($locusSection==1){
		push@loci,$line;
	}else{
		my@columns=split "\t", $line;
		shift@columns;
		foreach my$i (0..$#columns){
			my$allele1=substr $columns[$i], 0, 2;
			my$allele2=substr $columns[$i], 2, 2;
			if($allele1 eq "00"){
				$missingGenos{$loci[$i]}++;
			}elsif($allele1 ne $allele2){
				$hetGenos{$loci[$i]}++;
				$totalGenos{$loci[$i]}++;
			}else{
				$totalGenos{$loci[$i]}++;
			}
		}
	}
}

foreach my$key (sort keys %totalGenos){
	my$missingPerc;
	my$hetPerc;
	if(not exists $missingGenos{$key}){
		$missingGenos{$key}=0;
		$missingPerc=0
	}else{
		$missingPerc=$missingGenos{$key}/$totalGenos{$key};
	}
	if(not exists $hetGenos{$key}){
		$hetGenos{$key}=0;
		$hetPerc=0;
	}else{
		$hetPerc=$hetGenos{$key}/$totalGenos{$key};
	}
	print "$pop\t$key\t$totalGenos{$key}\t$missingGenos{$key}\t$hetGenos{$key}\t$missingPerc\t$hetPerc\n";
}
