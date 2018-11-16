#!usr/bin/perl

# Xin Qiao, 04/26/2017

use warnings;
use strict;

die "perl $0 collinearity.kaks\n" if $#ARGV < 0;

my ($sp) = ($ARGV[0] =~ /^(.{3})/);

my %sco;
my %eva;
my %siz;
my %chro;
my %ori;
open IN, "$ARGV[0]" or die "Cannot open $ARGV[0]!";
while(<IN>)
{
	chomp;
	if($_ =~ /^\#\# Alignment ([\d]+)\:/)
	{
		my $bkid="$1";
		my @a=split /\s/, $_;
		my @b=split /\=/, $a[3];
		$sco{$bkid}=$b[1];
		my @c=split /\=/, $a[4];
		$eva{$bkid}=$c[1];
		my @d=split /\=/, $a[5];
		$siz{$bkid}=$d[1];
		$chro{$bkid}=$a[6];
		$ori{$bkid}=$a[7];
	}
}
close IN;

my %bks;
open IN, "$ARGV[0]" or die "Cannot open $ARGV[0]!";
while(<IN>)
{
	chomp;
	if($_ =~ /^\#/){next;}
	$_ =~ s/^\s+//g;
	my @a=split /\t/, $_;
	if($a[5] eq "NA"){next;}
	if($a[5] > 5){next;}#remove the Ks values > 5.0
	my @b=split /\-/, $a[0];
	if(exists $bks{$b[0]}{$a[5]})
	{
		next;
	}
	else
	{
		$bks{$b[0]}{$a[5]}='A';
	}
}
close IN;

my $i=0;
my $sumks=0;
my $aveks;
my %baks;
foreach my $key1 (keys %bks)
{
	foreach my $key2 (keys %{$bks{$key1}})
	{
		$i++;
		$sumks += $key2;
	}
	$aveks = $sumks / $i;
	$baks{$key1}=$aveks;
	$i=0;
	$sumks=0;
}

open OUT, ">$sp.synteny.blocks.ks.info";
print OUT "Blocks ID\tLocation\tBlock Size\tAverage Ks\te-value\tScore\tOrientation\n";
foreach my $key1 (keys %bks)
{
	if(!exists $chro{$key1})
	{
		print "$key1\n";
	}
	else
	{
		my @a=split /\&/, $chro{$key1};
		#if($a[0] eq $a[1]){next;}#herein, we removed those intra-chromsomal synteny blocks
		#else
		#{
			print OUT "Alignment$key1\t$chro{$key1}\t$siz{$key1}\t$baks{$key1}\t$eva{$key1}\t$sco{$key1}\t$ori{$key1}\n";
		#}
	}
}
close OUT;

__END__