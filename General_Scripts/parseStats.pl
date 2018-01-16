#!/usr/bin/perl -w
use strict;

#Hash of single and duplicate marker models, key is model, value is 1 for single locus or 2 for duplicate.
my%dupModels=(	"AA_xx" => 1,
				"AB" => 1,
				"AA_AB" => 2,
				"AA_BB" => 2,
				"AB_AB" => 2,
				"AA_AB_AB" => 2,
				"AA_BB_AB" => 2,
				"AB_AB_AB" => 2,
				"AA_BC" => 2,
				"AB_AC" => 2,
				"AB_CD" => 2,
				"AA_BB_AC" => 2,
				"AA_AB_AC" => 2,
				"AA_BB_CC" => 2);


while(my$line=<>){
	chomp$line;
	my($tag,$mappable_loci,$epsilon,$bestModel,$bestLikelihood,$secondModel,$secondLikelihood,$segTest_x1,$segTest_x2)=split "\t", $line, 9;
	if($mappable_loci>=1){
		print "$tag\t$bestModel\t$dupModels{$bestModel}\n";
	}
}