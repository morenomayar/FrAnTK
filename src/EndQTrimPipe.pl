
use Getopt::Long;
my %opts=();
GetOptions (\%opts,'i=s', 'trim=i', 'bed=s');

open(BAM,"samtools view -L $opts{bed} -hF 4 $opts{i} |");
open(OUT, "| samtools view -bS -");

$trim=$opts{trim};
while(<BAM>){
if($_ =~ /^\@/){
        print OUT $_;
}else{
        $linec=$_;
        chomp($linec);
        @line=split(/\t/, $linec);
        $q=$line[10];
        @qual=split(//, $q);
        $length=$qual;
        for($i=0; $i<$trim; $i++){
                $qual[$i]="\"";
                $qual[$length-($i+1)]="\"";
        }
        $q=join("", @qual);
        $line[10]=$q;
        $samline=join("\t", @line);
        print OUT "$samline\n";
}
}
close(BAM);
close(OUT);

