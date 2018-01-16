#parsing genepop files to get obs hets and Fis 


#this script only works if you edit the DIV file first. Replace Locus with {space}Locus and replace the # of spaces until there is only one separating each of the columns with the hets and Fis. 

#usage

#perl genepop....pl genepop_al_freq_file.txt.div popnames.txt > out.txt
#argument 1: a .div file from genepop
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
      $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; 
      $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; $_=<IN>; 	  #skip 25 lines
       
      #if the data for this snp is not missing meaning it has ------ on the line after pop \t alleles\t genes
      if($_ =~ m/-------------------------------------------/)
      {
      #$_=<IN>; #skip a line
      @Qintra=();
      @Qinter=();
      @Fis=();
      #this may not work well for snps that are monomorphic in all pops but we don't want them anyway
      for($i=1; $i<=$num_pops; $i=$i+1) #iterates through each pop for each snp snp grabbing hets and fis
      {
        $freq_line=<IN>;
        @split_freq_line=split(" " , $freq_line); #split by spaces 
        $each_quintra="$split_freq_line[1]\t"; #grabbing the freq of al1
        $each_quinter="$split_freq_line[2]\t"; #grabbing the counts
        $each_fis="$split_freq_line[3]\t"; #grabbing the counts
        push(@Qintra,$each_quintra); 
        push(@Qinter,$each_quinter);
        push(@Fis,$each_fis);
		
      }
      
      print("@Qintra");
      print("@Qinter");
      print("@Fis\n");
           
     }
     else #if there is missing data just print snpid with nothing else
      {
        print("$title_line[1]\n");
        redo;
      }
     }
}
close IN;







