#Subset a specified set of SNPs and output appended genepop file

#arg1 SNPlist each line the full name of a SNP from the header in the genepop file
#arg2 genepop file

#usage
#perl subset_genepop_by_SNPs.pl SNP_names.txt 96_snps_2_1_genepop.txt > new_genepop_file.txt
#need to change this to match your SNP name prefix i.e. Ots, MY, One, RAD

print("Title Line:\n");

open(SNPs, "<$ARGV[0]") or die "Error!!! reading $ARGV[0] for reading";
@SNP_names=();
while(<SNPs>)
{
  chomp($_);
  $SNP="$_";
  push(@SNP_names, $SNP);
}
close SNPs;

#print("@SNP_names");


open(IN, "<$ARGV[1]") or die "Error!!! reading $ARGV[1] for reading"; #opens genepop file
$z=0;
@SNP_index=();
while(<IN>)
 {
    chomp($_);
    $line=$_;
    foreach $SNP (@SNP_names)
    {
      if($line eq $SNP) 
      {
      print("$line\n");     
      push(@SNP_index,$z);
      }      
    }
    if($line !~m/Title|Pop|pop|,|Stacks/){$z++;}
  }   
  close IN;  

open(IN, "<$ARGV[1]") or die "Error!!! reading $ARGV[1] for reading"; #opens genepop file
while(<IN>)
  {
  chomp($_);
  $line=$_;
    if($line=~m/Pop|pop/){print("Pop\n");}
    if($line =~ m/,/)
    {
    $ind=$line;
    @split_ind_line=split(" " , $ind);
    print("$split_ind_line[0]\t");
      foreach $index (@SNP_index){print("$split_ind_line[$index+1]\t");}
    print("\n");
    }
  }
  close IN;
  
  







