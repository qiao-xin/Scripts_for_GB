#!/usr/bin/perl

# Leiting Li, Oct 15, 2014
# Xin Qiao, Oct 05, 2018
# Xin Qiao, Sep 25, 2019

use warnings;
use strict;
use Bio::Perl;
use Bio::DB::Fasta;
use Getopt::Long;
use File::Basename;

sub usage{
	die "Usage: ", basename($0)," -d <CDS.fasta> -g <Gene_pairs> -o <outfile_name>\n";
}


die &usage if (@ARGV == 0 and -t STDIN);

my $db_file;
my $gene_file;
my $out_file;

my $result = GetOptions("d=s" => \$db_file, 
                        "g=s" => \$gene_file,
                        "o=s" => \$out_file);

&usage unless $db_file and  $gene_file and $out_file;
#print "DB file: $db_file\nGene file: $gene_file\nOut dir: $out_dir\n";

my $db = Bio::DB::Fasta->new($db_file);
open my $in_fh, $gene_file or die "$!";
my $out_dir = "tmp_dir.".$$.".".time;
die "Making temporary directory failed!\n" unless mkdir $out_dir;

my @contents = <$in_fh>;
@contents = grep{!/^\s*$/}@contents;
@contents = grep{!/^Duplicate/}@contents;
@contents = map{s/[\r\n]//g;$_}@contents;
my $amount = @contents;
my $num_len = length($amount);


for(my $count = 0; $count < @contents; $count++){
	my $out_pre = sprintf "$out_dir/%0${num_len}s", $count + 1;
	my $out = Bio::SeqIO->new(
			-file => "> $out_pre.fasta", 
			-format => 'fasta');
	my $str = $contents[$count];
	#my @genes = split /[,\s]+/, $str;
	my @genes = split /\t/, $str;
	my @genes2=();
	push(@genes2,$genes[0]);
	push(@genes2,$genes[2]);
	my $seq1 = $db->get_Seq_by_id($genes[0]);
	my $seq2 = $db->get_Seq_by_id($genes[2]);
	if(!$seq1 || !$seq2)
	{
		next;
	}
	for(@genes2){
		my $seq = $db->get_Seq_by_id($_);
		$out->write_seq($seq);
	}
	$out->close;
	system "perl fa_prepareDATA4_kaks_calculator.pl $out_pre.fasta > $out_pre.axt";
}

#my ($sp) = ($gene_file =~ /^(\S+?)\./);
#system "cat $out_dir/*.axt > $out_file.axt";#modified Sep 25 2019
system "find $out_dir/ -name \"*.axt\" | xargs cat > $out_file.axt.raw";#modified Sep 25 2019, added at 26 Oct 2020
system "mv $out_file.axt.raw $out_file.axt";#modified Sep 25 2019, added at 26 Oct 2020

warn "CAUTION: deleting temporary files in $out_dir failed!\n" unless unlink glob "$out_dir/*";
warn "CAUTION: removing temporary directory $out_dir failed!\n" unless rmdir $out_dir;

##################
#print "$out_file\n";

system "KaKs_Calculator -i $out_file.axt -m GMYN -o $out_file.KKC.format > $out_file.KKC.logfile";#modified Sep 25 2019

open IN, "$out_file.KKC.format";#modified Sep 25 2019
open OUT, ">$out_file";#modified Sep 25 2019
print OUT "Duplicate 1\tDuplicate 2\tKa\tKs\tKa/Ks\tP-Value\n";
while(<IN>)
{
	chomp;
	my @a=split /\t/, $_;
	my @b;
	@b=split /-/, $a[0];
	print OUT "$b[0]\t$b[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\n";
}
close OUT;
close IN;

system "sed -i '1i\\Gene pairs\tMethod\tKa\tKs\tKa/Ks\tP-Value (Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites (0:2:4)\tSubstitutions\t \\
S-Substitutions\tN-Substitutions\tFold-S-Substitutions (0:2:4)\tFold-N-Substitutions (0:2:4)\tDivergence-Time\tSubstitution-Rate-Ratio (rTC:rAG:rTA:rCG:rTG:rCA/rCA)\t \\
GC(1:2:3)\tML-Score\tAICc\tAkaike-Weight\tModel\n' $out_file.KKC.format";

__END__
