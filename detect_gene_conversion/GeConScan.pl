#!usr/bin/perl

# Xin Qiao, 12/29/2016

use warnings;
use strict;

die "Usage: perl $0 <Species_A_abbr(target)> <Species_B_abbr(outgroup)> <order(f or r)> <gene_pairs> \nIf the target genome gene_id located before the outgroup genome gene_id, please use 'f'\n" if $#ARGV < 3;

use Bio::SeqIO;
use Bio::Seq;

sub format_cds;

############ Parsing ortholog pairs between A and B genome
my %orp;
open IN, $ARGV[0]."_$ARGV[1].collinearity" or die "Cannot open ".$ARGV[0]."_$ARGV[1].collinearity!";
while(<IN>)
{
	chomp;
	if($_ =~ /^\#/){next;}
	my @a=split /\t/, $_;
	$a[3] =~ s/\s+//g;
	if($ARGV[2] eq 'f')
	{
		if(exists $orp{$a[1]}{$a[2]})
		{
			next;
		}
		else
		{
			$orp{$a[1]}{$a[2]}=$a[3];
		}
	}
	elsif($ARGV[2] eq 'r')
	{
		if(exists $orp{$a[2]}{$a[1]})
		{
			next;
		}
		else
		{
			$orp{$a[2]}{$a[1]}=$a[3];
		}
	}
}
close IN;

my %eva;
my %temp;
foreach my $key1 (keys %orp)
{
	foreach my $key2 (keys %{$orp{$key1}})
	{
		if(!exists $eva{$key1})
		{
			$eva{$key1}=$key2;
			$temp{$key1}=${$orp{$key1}}{$key2};
		}
		else
		{
			if(${$orp{$key1}}{$key2} < $temp{$key1})
			{
				$eva{$key1}=$key2;
				$temp{$key1}=${$orp{$key1}}{$key2};
			}
		}
	}
}
#############searching the 'A1 B1 A2 B2' quartets#########################
system 'mkdir alignment';
system 'mkdir alignment/clustdna';
system 'mkdir alignment/clustprot';
system 'mkdir alignment/pamldna';
open IN, $ARGV[3] or die "Cannot open $ARGV[3]!\n";
open OUT, ">".$ARGV[0]."_$ARGV[1].homologous.quartets";
print OUT "Duplicate 1	Homolog 1	Duplicate 2	Homolog 2\n";
while(<IN>)
{
	chomp;
	if($_ =~ /^Duplicate/){next;}
	my @a=split /\t/, $_;
	if(exists $eva{$a[0]} && exists $eva{$a[2]})
	{
		if($eva{$a[0]} ne $eva{$a[2]})
		{
			print OUT "$ARGV[0].$a[0]\t$ARGV[1].$eva{$a[0]}\t$ARGV[0].$a[2]\t$ARGV[1].$eva{$a[2]}\n";
		}
	}
	else
	{
		next;
	}
}
close OUT;
close IN;
##################### preparing CDS sequences
format_cds($ARGV[0]);
format_cds($ARGV[1]);
sub format_cds()
{
	my $abr = $_[0];
	open OUT, ">$abr.cds.temp";
	my $readseq = Bio::SeqIO->new(-file=>"$abr.cds", -format=>'fasta');
	while(my $seq_obj = $readseq->next_seq)
	{
		my $display_name = $seq_obj->display_name;
		my $seq = $seq_obj->seq;
		print OUT ">$abr.$display_name\n$seq\n";
	}
	close OUT;
	my $out1 = Bio::SeqIO->new(-file => ">$abr.cds.formatted", -format =>'fasta');
	my $readseq2 = Bio::SeqIO->new(-file=>"$abr.cds.temp", -format=>'fasta');
	while(my $seq_obj = $readseq2->next_seq)
	{
		my $display_name = $seq_obj->display_name;
		#my $seq_length = $seq_obj->length;
		my $seq = $seq_obj->seq;
		$out1->write_seq($seq_obj);
	}
	system "rm -f $abr.cds.temp";
}
system 'mkdir cds';
system "perl split.fasta.pl $ARGV[0].cds.formatted cds";
system "perl split.fasta.pl $ARGV[1].cds.formatted cds";
system "rm -f $ARGV[0].cds.formatted $ARGV[1].cds.formatted";
##################### detecting gene conversion events
my $out = $ARGV[0]."_$ARGV[1].homologous.quartets";
system "perl calculate.K.quartet.genes.CV.BP.all.pl $out $ARGV[0] $ARGV[1]";
open IN, "$ARGV[0]_$ARGV[1].homologous.quartets.K.CV.BP.txt";
open OUT, ">$ARGV[0]_$ARGV[1].quartets.K.CV.BP.txt";
print OUT "Note:	".$ARGV[0]."-1 and ".$ARGV[0]."-2 show two ".$ARGV[0]." paralogous genes, ".$ARGV[1]."-1 and ".$ARGV[1]."-2 show two ".$ARGV[1]." paralogous genes\n";
print OUT "	Pa() shows nonsynonymous nucleotide substitution difference, Ps() shows synonymous nucleotide substitution difference\n";
print OUT "	'N' shows no gene converison found between the $ARGV[0] or $ARGV[1] paralogs\n";
print OUT "	'R' shows that the two $ARGV[0] paralogs have been affected by gene covnersion\n";
print OUT "	'S' shows that the two $ARGV[1] paralogs have been affected by gene covnersion\n";
my $t1="Homologous_quartets\tPa($ARGV[0]-1, $ARGV[0]-2)\tPs($ARGV[0]-1, $ARGV[0]-2)\tPa($ARGV[1]-1, $ARGV[1]-2)\t";
my $t2="Ps($ARGV[1]-1, $ARGV[1]-2)\tPa($ARGV[0]-1, $ARGV[1]-1)\tPs($ARGV[0]-1, $ARGV[1]-1)\tPa($ARGV[0]-2, $ARGV[1]-2)\t";
my $t3="Ps($ARGV[0]-2, $ARGV[1]-2)\tIs_conversion_between_$ARGV[0]-1_and_$ARGV[0]-2\tIs_conversion_between_$ARGV[1]-1_and_$ARGV[1]-2\t";
my $t4="Bootstrap_percentage_Ks_$ARGV[0]-1_$ARGV[0]-2\tBootstrap_percentage_Ka_$ARGV[0]-1_$ARGV[0]-2\t";
my $t5="Bootstrap_percentage_Ks_$ARGV[1]-1_$ARGV[1]-2\tBootstrap_percentage_Ka_$ARGV[1]-1_$ARGV[1]-2";
print OUT "$t1$t2$t3$t4$t5\n";
while(<IN>)
{
	if($_ =~ /^groupsize/)
	{
		next;
	}
	else
	{
		print OUT "$_";
	}
}
close OUT;
close IN;
system "rm -f $ARGV[0]_$ARGV[1].homologous.quartets.K.CV.BP.txt";
system 'rm -f yn00.ctl';
system 'rm -rf alignment';
system 'rm -rf cds';
system 'rm -f 2YN.dN 2YN.dS 2YN.t';
system 'rm -f rst rst1 rub';

__END__