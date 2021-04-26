
  # mpileup format:
  # chr pos             Ref     Depth1  Pileup1 Qual1   Depth2  Pileup2 Qual2
  # 22  32781587        N       1       ^FG     G       1       ^FG     G
  # 22  32781722        N       1       c       E       1       c       E
  # 22  32904927        N       2       t$t     HB      0       *       *
  # 22  42815298        N       0                       0       *       *
  # 22  42815298        N       1       c       H       0       ^~A$    J



sub parse_pileup {
  # Input:
  #   Pileup column (e.g: '^FG' 'c' 't$t');
  # Output:
  #   (number_of_A,number_of_C,number_of_G,number_of_T)
  my $PU = shift;
  $PU = uc($PU);
  my %alleles=();
  $alleles{"A"}=0;
  $alleles{"C"}=0;
  $alleles{"G"}=0;
  $alleles{"T"}=0;
  # A pileup line may look like this G^FG^Fg^FG^FG^FG+1GC-2NN*GA$G
  while($PU) {
    #Always get the first operation into $1 and the rest of the operations into $2
    if($PU =~ s/^(\w)//) {
      # E.g PU="G"
      # Whenever there is a match, add that match to the corresponding allele count
      $alleles{$1}++;
      next;
    } elsif($PU =~ s/^\^(.)//) {
      # E.g PU="^F"
      # When a read ends, there is a hat followed by the read mapping quality. Get rid of this info.
      next;
    } elsif($PU =~ s/^\$//) {
      # E.g PU='$'
      # A dollar means a read starts here, get rid of it.
      next;
    } elsif($PU =~ s/^\*//) {
      # E.g PU='*'
      # * represents a deleted base in a previous read - ignore that
      next;
    } elsif($PU =~ s/^[-+](\d+)(.*)//) {
      # E.g. -2NN
      # This one deals with indels. +2GG would mean that 2 Gs were inserted between this position and the ref in one read
      my $nInDel = $1;
      # Skip the appropriate number of positions in the remaining pileup line
      $PU =~ s/^\w{$nInDel}//;
    } else {
      die_bug("PU contains unmatched chars: $PU");
    }
  }
  # Return the array of allele counts
  return @alleles{"A","C","G","T"};
}


while(<STDIN>){
	$line=$_;
	chomp($line);
	@l=split(/\t/, $line);
	@AC=parse_pileup($l[4]);
	$AlleleString="A"x$AC[0]."C"x$AC[1]."G"x$AC[2]."T"x$AC[3];
	$Alength=length($AlleleString);
	@AS=split(//, $AlleleString);
	$randomallele= $AS[rand @AS];
	if($l[3]!=0 && $Alength!=0){
		print $l[0]."_".$l[1]."\t".$randomallele."\n";
	}
}



















