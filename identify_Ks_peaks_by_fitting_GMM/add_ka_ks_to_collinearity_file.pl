#!usr/bin/perl

# Xin Qiao, 10/11/2018

use warnings;
use strict;

die "perl $0 species_abbr\n" if $#ARGV < 0;

sub read_kaks;
sub add_kaks;

my %ka;
my %ks;
my %kaks;
read_kaks("$ARGV[0].wgd.kaks",$ARGV[0]);
add_kaks($ARGV[0]);

sub read_kaks()
{
	my $f = $_[0];
	open IN, $f or die "Cannot open $f!\n";
	while(<IN>)
	{
		chomp;
		my @a=split /\t/, $_;
		my $pr="$a[0]\t$a[1]";
		$ka{$pr}=$a[2];
		$ks{$pr}=$a[3];
		$kaks{$pr}=$a[4];
	}
	close IN;
}

sub add_kaks()
{
	my $sp = $_[0];
	open IN, "$sp.collinearity" or die "Cannot open $sp.collinearity!\n";
	open OUT, ">$sp.collinearity.kaks";
	while(<IN>)
	{
		chomp;
		if($_ =~ /^\#/)
		{
			print OUT "$_\n";
		}
		else
		{
			my @a=split /\t/, $_;
			my $pr="$a[1]\t$a[2]";
			if(!exists $ka{$pr} && !exists $ks{$pr})
			{
				print OUT "$_\tNA\tNA\tNA\n";
			}
			else
			{
				print OUT "$_\t$ka{$pr}\t$ks{$pr}\t$kaks{$pr}\n";
			}
		}
	}
	close OUT;
	close IN;
}

__END__