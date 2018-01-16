#/usr/bin/perl -w
use strict;

foreach my$i (1..6){
	my$command1="populations -b 17 -P stacks -M stacks/popMap.txt -W stacks/whitelist_".$i.".txt --genepop";
	my$command2="mkdir stacks/whitelist_".$i;
	my$command3="cp stacks/batch_17* stacks/whitelist_".$i;
	print "$command1\n";
	`$command1`;
	print "$command2\n";
	`$command2`;
	print "$command3\n";
	`$command3`;
}
