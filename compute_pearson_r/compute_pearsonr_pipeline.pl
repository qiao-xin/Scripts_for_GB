#!usr/bin/perl

# Xin Qiao, 01/23/2018
# Xin Qiao, 01/16/2019

use warnings;
use strict;

# Before running this pipeline, you should run Kallisto software (https://github.com/pachterlab/kallisto) to quantify expression abundances for each gene.
# This script was used to integrate the expression profiles of different development stages or conditions for each gene, 
# and then to compute Pearson correlation coefficients (r) between expression profiles of the two gene copies.
# The python script 'pearsonr.py' is required for 

#die "Usage: perl $0 <>\n" if $#ARGV < 0;

sub log10;

my %cutoff;
$cutoff{"Ath"}=0.02;#just for example

#my %dirn;
my %expr;
print "Loading TPM values from Kallisto results...\n";
print "================================================\n\n";
my %id=();
my %anc=();
my %anc2=();
system "ls ./example_data -lR|grep \"^d\"|wc -l >dir.num";
my $dirn=0;
open IN, "dir.num";
while(<IN>)
{
	chomp;
	#$dirn{$key}=$_;
	$dirn=$_;
}
close IN;
system "rm dir.num";
system "cp ./example_data/output1/abundance.tsv ./";
open IN, "abundance.tsv" or die "Cannot open abundance.tsv!\n";
while(<IN>)
{
	chomp;
	if($_ =~ /^target\_id/){next;}
	my @a=split /\t/, $_;
	$id{$a[0]}='A';
}
system 'rm abundance.tsv';
for(my $i=1; $i<=$dirn; $i++)
{
	#print $i, "\n";
	system "cp ./example_data/output$i/abundance.tsv ./";
	open IN, "abundance.tsv" or die "Cannot open abundance.tsv!\n";
	while(<IN>)
	{
		chomp;
		if($_ =~ /^target\_id/){next;}
		my @a=split /\t/, $_;
		my $info="";
		$info="$i$a[0]";
		if($a[4] eq "-nan")
		{
			$anc2{$info}=0;
			$anc{$info}=log10(0.0001);
		}
		elsif($a[4] eq "nan")
		{
			$anc2{$info}=0;
			$anc{$info}=log10(0.0001);
		}
		elsif($a[4] == 0)
		{
			$anc2{$info}=$a[4];
			$anc{$info}=log10(0.0001);
		}
		else
		{
			$anc2{$info}=$a[4];
			$anc{$info}=log10($a[4]);
		}
	}
	system 'rm abundance.tsv';
}
open OUT, ">alltpm";
open OUT2, ">alltpm.log10";
my $out="GeneID";
for(my $i=1; $i<=$dirn; $i++)
{
	$out .= "\t$i";
}
print OUT "$out\n";
print OUT2 "$out\n";
foreach my $key (keys %id)
{
	print OUT "$key\t";
	print OUT2 "$key\t";
	for(my $i=1; $i<=$dirn; $i++)
	{
		my $b="$i$key";
		if($i < $dirn)
		{
			print OUT "$anc2{$b}\t";
			print OUT2 "$anc{$b}\t";
		}
		elsif($i = $dirn)
		{
			print OUT "$anc2{$b}\n";
			print OUT2 "$anc{$b}\n";
		}
	}
}
close OUT;
close OUT2;
	
################## Reading the sum of expression vaules for each gene in different conditions
#my %expr=();
my $sum=0;
open IN, "alltpm" or die "Cannot open alltpm!\n";
while(<IN>)
{
	chomp;
	if($_ =~ /^GeneID/){next;}
	my @a=split /\t/, $_;
	for(my $i=1; $i<=$dirn; $i++)
	{
		$sum += $a[$i];
	}
	$expr{$a[0]}=$sum;
	$sum=0;
}
close IN;
	
################## Generate 10,000 randomly selected gene pairs
my $g=0;
my $gpair="";
open OUT, ">random.pairs";
my $n=0;
foreach my $key (keys %id)
{
	$n++;
	$g++;
	if($g ==1)
	{
		$gpair .= "$key\t";
	}
	elsif($g == 2)
	{
		$gpair .= "$key";
		print OUT "$gpair\n";
		$g=0;
		$gpair="";
		#last;
	}
	if($n == 20000)
	{
		last;
	}
}
close OUT;
#%expr=();

my @files = glob "*.pairs";
foreach my $pairs (@files)
{
	print "Parsing $pairs ...\n";
	print "---------------------------------------\n\n";
	if($pairs =~ /random/)
	{
		open IN, "random.pairs" or die "Cannot open random.pairs!\n";
		open OUT, ">random.pairs.filtered";
		while(<IN>)
		{
			chomp;
			my @a=split /\t/, $_;
			if($expr{$a[0]} >= $cutoff{$sp} && $expr{$a[1]} >= $cutoff{$sp})
			{
				print OUT "$_\n";
			}
			else
			{
				next;
			}
		}
		close OUT;
		close IN;
		system "rm -f random.pairs";
		print "   Estimating the pearson r...\n";
		print "   -----------------------------------------\n\n";
		system "python pearsonr.py alltpm.log10 random.pairs.filtered $pairs.pr";
		system "rm -f random.pairs.filtered";
	}
	else
	{
		system "cat $pairs | perl -pe \'s/^Duplicate.*\\n//g;\' | perl -pe \'s/^(\\S+\\t)\\S+\\t(\\S+)\\t\\S+\\t\\S+\\n/\$1\$2\\n/g;\' | perl -pe \'s/^(\\S+\\t)\\S+\\t(\\S+)\\t\\S+\\t\\n/\$1\$2\\n/g;\' >$pairs.raw";
		open IN, "$pairs.raw" or die "Cannot open $pairs.raw!\n";
		open OUT, ">$pairs.new";
		while(<IN>)
		{
			chomp;
			my @a=split /\t/, $_;
			if(!exists $expr{$a[0]} || !exists $expr{$a[1]})
			{
				next;
			}
			else
			{
			if($expr{$a[0]} >= 0.01 && $expr{$a[1]} >= 0.01)#selected gene pairs in which both copies are expressed
			                              #in at least one tissue (TPM >= cutoff, herein, 0.01 was used for example)
			{
				print OUT "$_\n";
			}
			else
			{
				next;
			}
			}
		}
		close OUT;
		close IN;
		system "rm -f $pairs.raw";
		print "   Estimating the pearson r...\n";
		print "   -----------------------------------------\n\n";
		system "python pearsonr.py alltpm.log10 $pairs.new $pairs.pr";
		system "rm -f $pairs.new";
	}
}

sub log10()
{
	my $n = shift;
	return log($n)/log(10);
}

__END__