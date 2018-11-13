#!/usr/bin/perl

# Leiting Li, Apr 25, 2013
# Leiting Li, May 7, 2013
# Leiting Li, Sep 15, 2014
# Xin Qiao, Oct 5, 2018

# Start
use warnings;
use strict;
use Bio::Perl;
use Getopt::Long;
use File::Basename;

my $model = "P";
my $outfile = '';
my $strict_CDS;
my $verbose;
GetOptions ("model=s" => \$model, 
			"out=s" => \$outfile,
			"strict" => \$strict_CDS,
            "verbose" => \$verbose);

if( @ARGV == 0 and -t STDIN ){ &usage; }
my $file_cds;
if(@ARGV){
	$file_cds = shift @ARGV;
}else{
	$file_cds = "-";
}
my $debug = 0;

#
# Usage
#
# -m, --model=Model    
#    P for pairwise, S for sequence, default:P

sub usage{
    die "\n",basename($0)," [OPTINOS] <CDS>

OPTIONS
  -o, --out=FILE    
      output file name, defaul: STDOUT

  -s, --strict        
      strict CDS, enable this option will check  start 
      codon, stop codon, and the integrity of codons.

  -v, --verbose
      print verbose information

  Note: CDS should be in FASTA format
  CAUTION: sequence name should not contain signals \\, |, /, etc \n\n";
}

&main;

sub main{
	my @all_cds = read_all_sequences($file_cds,'fasta');
	&check_cds(\@all_cds);
	@all_cds = &mask_stop_codon(\@all_cds);
	die unless @all_cds;

	my $tmp_dir = &make_tmp_dir;
	my $out_fh;
	if($outfile){
		if($outfile eq $file_cds){die "CAUTION: Output file name \"$outfile\" is the same as CDS file name!\n"}
		open $out_fh, "> $outfile" or die "$!";
	}else{
		$out_fh = \*STDOUT;
	}

	&cal_kaks(\@all_cds, $tmp_dir, $out_fh);
}

sub make_tmp_dir{
	my $tmp_dir =  "tmpdir-axt-converter-$$-".time;
	if(-e $tmp_dir){
		unlink (glob "$tmp_dir/*");
	}elsif(-w "."){
		mkdir $tmp_dir;
	}else{
		$tmp_dir = "/tmp/$tmp_dir";
		mkdir $tmp_dir;
	}
	return $tmp_dir;
}

sub check_cds{
	my $seqArray = shift @_;
	for(my $i = 0; $i <= $#{$seqArray}; $i++){
		$_ = $seqArray->[$i];
		my $length = $_->length;
		my $start_codon = $_->subseq(1,3);
		my $end_codon = $_->subseq($length - 2, $length);

		die unless $length and $start_codon and $end_codon;

		warn "++ Checking gene ",$i+1,": ", $_->display_id, " ... \n" if $verbose;
		if($length % 3 == 0){
			warn "++ \tNuleotide length is all right!\n" if $verbose;
		}else{
			die $_->display_id, ": Nucleotie length ERROR! $length % 3 = ".($length % 3)."!\n";
		}

		if($strict_CDS){
			if($start_codon =~ /^ATG$/i and $end_codon =~ /TAA|TAG|TGA/i){
				warn "++ \tLength, start codon, end codon OK!\n" if $verbose; 
			}else{
				die "Coding sequence ERROR, Length = $length, Start codon = $start_codon, End codon = $end_codon\n";
			}

			if(&with_stop_codon($_)){
				die $_->display_id, ": The coding sequence has stop codon not in the end!\n";
			}else{
				warn "++ \tNo stop codon except the end codon!\n" if $verbose;
			}
		}
	}
}

sub get_codon{
	my $seqstr = shift;
	my @nts = split //, $seqstr;
	my @codons;
	for(my $i = 0; $i <= $#nts; $i += 3){
		my $codon = join('', @nts[$i, $i+1, $i+2]);
		push @codons, $codon;
	}
	return @codons;
}

sub with_stop_codon{
	my $seqobj = shift @_;
	my @codons = &get_codon($seqobj->seq);
	for(my $i = 0; $i < $#codons; $i++){
		return 1 if $codons[$i] =~ /TAA|TAG|TGA/;
	}
	return 0;
}

sub mask_stop_codon{
	my $seqArray = shift @_;
	my @newArray;
	die "Empty array" unless @$seqArray;

	for (my $i = 0; $i <= $#{$seqArray}; $i++){
		my $seqobj = Bio::Seq->new(
			-display_id => $seqArray->[$i]->display_id,
			-seq => join('', map{s/TAA|TAG|TGA/NNN/i;$_}(&get_codon($seqArray->[$i]->seq))),
			-desc => $seqArray->[$i]->desc,
			-alphabet => $seqArray->[$i]->alphabet
			);
		push @newArray, $seqobj;
	}
	die "Empty array" unless @newArray;
	return @newArray;
}

sub run_cmd{
	warn "+ @_\n" if $verbose;
	system @_;
}

sub cal_kaks{
	my $seqArray = shift @_;
	my $tmp_dir = shift @_;
	my $out_fh = shift @_;
	##########

	my $num_of_seq = @{$seqArray};
	my $num_of_combination = $num_of_seq * ($num_of_seq - 1);
	my $num_length = length($num_of_combination);

	my $c = 0;
	for(my $i = 0; $i < $#{$seqArray}; $i++){
		for (my $j = $i + 1; $j <= $#{$seqArray}; $j++){
			$c++;
			my $prefix = sprintf "$tmp_dir/c%0${num_length}d", $c;

			warn "++ :$c Proccessing sequences ", $seqArray->[$i]->display_id, " and ", $seqArray->[$j]->display_id, " ... \n" if $verbose;
			write_sequence("> $prefix.cds.fasta", 'fasta', $seqArray->[$i], $seqArray->[$j]);
			write_sequence("> $prefix.pep.fasta", 'fasta', $seqArray->[$i]->translate(), $seqArray->[$j]->translate());
			#&run_cmd("clustalw2 -INFILE=$prefix.pep.fasta -PWMATRIX=BLOSUM -TYPE=PROTEIN -OUTPUT=FASTA -OUTFILE=$prefix.pep.aln > $prefix.clustalw2.out");
			&run_cmd("mafft --localpair --maxiterate 1000 --quiet $prefix.pep.fasta >$prefix.pep.aln");
			&run_cmd("pal2nal.pl $prefix.pep.aln $prefix.cds.fasta -output fasta > $prefix.coding.aln");
			&run_cmd("perl parseFastaIntoAXT.pl $prefix.coding.aln > $prefix.parseFastaIntoAXT.out");
			#system "KaKs_Calculator -i $prefix.coding.aln.axt -m YN -o $prefix.coding.aln.kaks > $prefix.KaKs_Calculator.out";
			die "ERROR -- $prefix.pep.aln -- File size is zero!!!" unless -s "$prefix.coding.aln";

			open my $fh, "$prefix.coding.aln.axt" or die;
			while(my $line = <$fh>){
				print $out_fh $line;
			}
			close $fh;

			unlink glob "$prefix*";
		}
	}
	rmdir $tmp_dir;
}

__END__
