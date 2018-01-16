#parsing genepop files to get allele frequency and allele counts (2N)

#this script should be generalized and can probably be used for everything

#usage

#perl genepop....pl genepop_al_freq_file.txt.inf popnames.txt > out.txt
#argument 1: a .inf file from genepop
#argument 2: a file with your popnames in order, one on each line, with no spaces, see example file

print("Locus\t");

open(POPS, "<$ARGV[1]") or die "Error!!! reading $ARGV[1] for reading";
#makes a header with pop1 \t pop2\t pop1\t pop2\n 
@pop_names=();
while(<POPS>)
{
  chomp($_);
  $pop="$_\t";
  push(@pop_names, $pop);
}
close POPS;

print("@pop_names");
print("@pop_names\n");

$num_pops=@pop_names;



open(IN, "<$ARGV[0]") or die "Error!!! reading $ARGV[0] for reading"; #opens genepop file
while(<IN>)
 {
    chomp($_);
    $line=$_;
    if($line =~ m/^[:space:]*.Locus/) #looking for tables where lines start with Locus (towards end of file)
    {
      @title_line=split(/[:blank:]/ , $line); 
      chomp($title_line[1]);
      print("$title_line[1]\t"); #this is just grabbing the SNP id
      $_=<IN>; $_=<IN>; $_=<IN>; #skip 3 lines, this shouldn't have to be changed
       
      #if the data for this snp is not missing meaning it has ------ on the line after pop \t alleles\t genes
      if($_ =~ m/-----/) #ex of missing data is snp 427 in WAK snp discov is basically Locus: \n-----\n No data 
      {
      $_=<IN>; #skip a line
      @al_freqs=();
      @al_counts=();
      #this may not work well for snps that are monomorphic in all pops but we don't want them anyway
      for($i=1; $i<=$num_pops; $i=$i+1) #iterates through each pop for each snp snp grabbing freqs and counts
      {
        $freq_line=<IN>;
        @split_freq_line=split(" " , $freq_line); #split by spaces 
        $each_freq="$split_freq_line[1]\t"; #grabbing the freq of al1
        $each_count="$split_freq_line[3]\t"; #grabbing the counts
        push(@al_freqs,$each_freq); 
        push(@al_counts,$each_count);
        
      }
      
      print("@al_freqs");
      print("@al_counts\n");
           
     }
     else #if there is missing data just print snpid with nothing else
      {
        print("$title_line[1]\n");
        redo;
      }
     }
}
close IN;







