use strict;
use Bio::SeqIO;

#print " ========= Usage ======= \nperl split.fasta.pl inputfile outputdir\n=========================\n";

my $fastafile = $ARGV[0];
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $fastafile);

while(my $seqobj = $is -> next_seq())
{
   my $outfile = $ARGV[1]."/".$seqobj->display_id().".fasta";
   my $os = Bio::SeqIO -> new(-format => 'fasta', -file => ">".$outfile); 
   $os -> write_seq($seqobj);
}


