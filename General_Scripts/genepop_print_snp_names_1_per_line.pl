open(IN1, "<$ARGV[0]") or die "Error!!! reading $ARGV[0] for reading";

$z=0;
while(<IN1>)
{
  $line=$_;
  chomp($line);
  if($line=~m/Stacks/){print("Title Line:\n");}
  if($z eq 1)
  {
  @SNP_names=split(",",$line);
  foreach $SNP (@SNP_names){print("$SNP\n");}
  }
  elsif($z>1){print("$line\n");}
  
  $z++;
}
close IN1;