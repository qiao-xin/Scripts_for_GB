#!usr/bin/perl

# Developed by Dr. Xiyin Wang
# Revised by Xin Qiao, 12/29/2016

use warnings;
use strict;

use Bio::SeqIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Seq::EncodedSeq;
use Bio::AlignIO;

sub main;

my $BOOTNUM = 1000;
my $input = $ARGV[0];
open(IN, $input) or die "cannot open block files $input due to $!.\n";
my @arr = split(/\//, $input);
my $output = $arr[$#arr].".K.CV.BP.txt";
open(CV, ">".$output) or die "cannot open outfile $output due to $!.\n";
my $lineno = 0;
while(<IN>)
{
   if($lineno eq 0){$lineno++; next;}
   #$_ =~ s/[\n\r]//g;
   chomp;
   my @temarr = split(/\t/, $_);

   my @osid1 = split(/[;\/\|]/, $temarr[0]);
   for(my $i1=0; $i1<=$#osid1; $i1++)
   { 
      my @osid2 = split(/[;\/\|]/, $temarr[2]);
      for(my $i2=0; $i2<=$#osid2; $i2++)
      {
         my @sbid1 = split(/[;\/\|]/, $temarr[1]);
         for(my $i3=0; $i3<=$#sbid1; $i3++)
         {
            my @sbid2 = split(/[;\/\|]/, $temarr[3]);
            for(my $i4=0; $i4<=$#sbid2; $i4++)
            {
               print CV "groupsize ".$#osid1." ".$#osid2." ".$#sbid1." ".$#sbid2."\n";
               my $quartetid = $osid1[$i1]."_".$osid2[$i2]."_".$sbid1[$i3]."_".$sbid2[$i4];
               main($quartetid);
            }
         }
      }
   }
}

sub main()
{
   my $quartetid = $_[0];
   my @temarr = split(/_/, $quartetid);
   my ($osid1, $osid2, $sbid1, $sbid2) = @temarr;

   my @dnaobj1 = ();
   my %dna_hash1;
   my $protos1 = Bio::SeqIO -> new(-file=> ">prot1.fasta", -format=>"fasta");

   my $seqno = 0;
   my $isallfileexist = 1;

   for(my $i=0; $i<=$#temarr; $i++)
   {
      if($temarr[$i] !~ /^\w/){next;}
      if($temarr[$i] =~ /^$ARGV[1]/)
      {
         my $fasta = "cds/".$temarr[$i].".fasta";

         if(!(-e $fasta)){$isallfileexist = 0;}

         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);
         $dnaobj1[$seqno] = $fastaio -> next_seq;
         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
      }
      elsif($temarr[$i] =~ /^$ARGV[2]/)
      {
         my $fasta = "cds/".$temarr[$i].".fasta";

         if(!(-e $fasta)){$isallfileexist = 0;}

         my $fastaio = new Bio::SeqIO(-format => "FASTA", -file => $fasta);
         $dnaobj1[$seqno] = $fastaio -> next_seq;
         $dnaobj1[$seqno] -> display_id($dnaobj1[$seqno] -> display_id());
         $dna_hash1{$dnaobj1[$seqno]->display_id()} = $dnaobj1[$seqno];
         $protos1 -> write_seq($dnaobj1[$seqno] -> translate());
      }
      $seqno ++;
   }

   if($isallfileexist eq 0){next;}

###### without outgroup alignment
   system("clustalw -infile=prot1.fasta");

   my $is_prot_aln1 = Bio::AlignIO -> new(-file=>"prot1.aln", -format=>"CLUSTALW");
   my $prot_aln1 = $is_prot_aln1 -> next_aln();

   #### length of alignment
   my $len  = $prot_aln1 -> length;

   my @id;
   my @protseq;
   my $seq_no = $prot_aln1 -> num_sequences;

   #### read the id and sequence
   for(my $i=0; $i<$seq_no; $i++)
   {
       $id[$i] = $prot_aln1 -> get_seq_by_pos($i+1) -> display_id();
       my @tmpseq = split(//, $prot_aln1 -> get_seq_by_pos($i+1) -> seq());
       $protseq[$i] = \@tmpseq;
   }

   #### minimum aa identity among genes
   my $mini_identity = 1;
   my $max_gap_level = 0;
   my $max_gap_num   = 0;
   for(my $k=0; $k<$seq_no; $k++)
   {
      for(my $j=$k+1; $j<$seq_no; $j++)
      {
          my $gapnum = 0;
          my $idennum= 0;
          for(my $i=0; $i<$len; $i++)
          {
             if($protseq[$k][$i] eq "-" || $protseq[$j][$i] eq "-"){$gapnum ++;}
             if($protseq[$k][$i] eq $protseq[$j][$i]){$idennum ++;}
          }
          
          my $identity;
          if($len - $gapnum ne 0)
          {
            $identity = $idennum/($len-$gapnum);
          }
          else
          {
            $identity = 0;
          }
          my $gap_level = $gapnum/$len;
          if($identity < $mini_identity){$mini_identity = $identity;}
          if($gap_level > $max_gap_level){$max_gap_level = $gap_level;}
          if($gapnum > $max_gap_num){$max_gap_num = $gapnum;}
      }
   }

   #### a restriction to sequence similarity     $max_gap_level > 0.20 ||
   if($mini_identity < 0.40 ||  $len - $max_gap_num < 50){next;}

   system("rm prot1.fasta prot1.aln prot1.dnd");
   my $os_prot_aln = Bio::AlignIO -> new(-file=>">alignment/clustprot/".$quartetid.".prot.aln", -format=>"CLUSTALW");
   $os_prot_aln -> write_aln($prot_aln1);

   my $dna_aln = &aa_to_dna_aln($prot_aln1, \%dna_hash1);
   my $os_dna_aln = Bio::AlignIO -> new(-file=>">alignment/clustdna/".$quartetid.".dna.aln", -format=>"CLUSTALW");
   $os_dna_aln -> write_aln($dna_aln);

   my $cdsalnfile = "alignment/pamldna/".$quartetid.".NUC";
   write_aln_paml($cdsalnfile, $dna_aln);
   my $seqfile = "alignment/pamldna/".$quartetid.".NUC";
   my $outfile = "KV.yn00";
   set_ctl_parameters("seqfile", $seqfile,
                      "outfile", $outfile);
   system("yn00");

   my $os_para = $osid1."-".$osid2;
   my $sb_para = $sbid1."-".$sbid2;
   my $os_sb_ortho1 = $osid1."-".$sbid1;
   my $os_sb_ortho2 = $osid2."-".$sbid2;

   my (%KaKs) = read_KYang($os_para, $sb_para, $os_sb_ortho1, $os_sb_ortho2, "KV.yn00");
   unlink($seqfile, $outfile);

   my (%Ka, %Ks);
   foreach my $pair(keys(%KaKs))
   {
      my @pairvalue = split(/_/, $KaKs{$pair});
      $Ka{$pair} = $pairvalue[0];
      $Ks{$pair} = $pairvalue[1];
   }

   print  CV $quartetid."\t";
   print  CV $Ka{$os_para}."\t".$Ks{$os_para}."\t".$Ka{$sb_para}."\t".$Ks{$sb_para}."\t";
   print  CV $Ka{$os_sb_ortho1}."\t".$Ks{$os_sb_ortho1}."\t".$Ka{$os_sb_ortho2}."\t".$Ks{$os_sb_ortho2}."\t";

   if($Ks{$os_para} < 0 || $Ks{$sb_para} < 0 || $Ks{$os_sb_ortho1}  < 0 || $Ks{$os_sb_ortho2} < 0)
   {
      print CV "N\tN\n"; next;
   }

   my $is_conversion_os=0;
   my $is_conversion_sb=0;

   if($Ks{$os_para} <= $Ks{$os_sb_ortho1} && $Ks{$os_para} <= $Ks{$os_sb_ortho2})
   {
      $is_conversion_os = 1; print CV "R\t";
   }
   else
   {  print CV "N\t";
   }

   if($Ks{$sb_para} <= $Ks{$os_sb_ortho1} && $Ks{$sb_para} <= $Ks{$os_sb_ortho2})
   {
      $is_conversion_sb = 1; print CV "S\t";
   }
   else
   {  print CV "N\t";
   }

   my @bootprob;
   if($is_conversion_os +  $is_conversion_sb ne 0)
   {
      @bootprob = bootstrap_global($os_para, $sb_para, $os_sb_ortho1, $os_sb_ortho2, $dna_aln);
   }
   if($is_conversion_os eq 1 &&  $is_conversion_sb eq 1)
   {
     print CV "@bootprob[0..$#bootprob]\n";
   }
   elsif($is_conversion_os eq 1)
   {
     print CV "@bootprob[0..1] - -\n";
   }
   elsif($is_conversion_sb eq 1)
   {
     print CV "- - @bootprob[2..3]\n";
   }
   else
   {
     print CV "\n";
   }
   
   $lineno ++;
}
close($ARGV[0]);

sub read_KYang
{
   my $os_para = $_[0];
   my $sb_para = $_[1];
   my $os_sb_ortho1 = $_[2];
   my $os_sb_ortho2 = $_[3];
   my $file = $_[4];

   open(YN, $file) or die "cannot open yn00 output file \n"; #or die "cannot open $file due to $!\n";
   my (%KaKs);
   my $countline = -1;
   my @id = ();
   while(<YN>)
   {
      my $line = $_;
      if($line =~ /^Nei \& Gojobori/)
      {
         $countline = 0; next;
      }
      if($countline ne -1 && $countline <= 4 && ($line =~ /^$ARGV[1]/ || $line =~ /^$ARGV[2]/ ))
      {
        $countline ++;

         $line =~ s/[\n\r]//g;
         $line =~ s/[()]/ /g;
         $line =~ s/\s+/ /g;

#         print ".............".$line;
         my @tmparr = split(/ /, $line);
#         print "............>@tmparr[0..$#tmparr]\n";

         for(my $i=0; $i<=($countline - 1) * 3; $i = $i +3)
         {
            if($i eq 0){$id[$countline] = $tmparr[0];}
            else
            {
               my $id1 = $id[int($i/3)];
               my $id2 = $id[$countline];
               #print $id1." ".$id2." ".$tmparr[$i-1]." ".$tmparr[$i]."\n";
               $KaKs{$id1."-".$id2} = $tmparr[$i-1]."_".$tmparr[$i];
               $KaKs{$id2."-".$id1} = $tmparr[$i-1]."_".$tmparr[$i];
            }
         }
      }
      if($countline > 4){last;}
   }

   #print  "....".$os_para." ".$KaKs{$os_para}."\t".$sb_para." ".$KaKs{$sb_para}."\t";
   #print  "....".$os_sb_ortho1." ".$KaKs{$os_sb_ortho1}."\t".$os_sb_ortho2." ".$KaKs{$os_sb_ortho2}."\n";

   return ($os_para, $KaKs{$os_para},
           $sb_para, $KaKs{$sb_para},
           $os_sb_ortho1, $KaKs{$os_sb_ortho1},
           $os_sb_ortho2, $KaKs{$os_sb_ortho2}
   );
}

sub bootstrap_global()
{
   print "in bootstraping ................\n";
   my $os_para = $_[0];
   my $sb_para = $_[1];
   my $os_sb_ortho1 = $_[2];
   my $os_sb_ortho2 = $_[3];
   my $dna_aln = $_[4];

#   $os_check -> write_aln($dna_aln);

   my $dna_aln_rand = $dna_aln;

   my $align_len = $dna_aln -> length();
   my %align_seq;
   my @seq;
   my %align_seq_rand;
   my @tmpid;
   my $i = 0;
   foreach my $seq_obj ($dna_aln -> each_seq())
   {
      $align_seq{$seq_obj -> display_id()} = $seq_obj -> seq();
      my @arr = split(//, $seq_obj -> seq());
      $seq[$i] = \@arr;
#      print "\n############\n".$seq_obj -> seq()."\n";
      $align_seq_rand{$seq_obj -> display_id()} = "";
      $tmpid[$i] = $seq_obj -> display_id();
      $i = $i + 1;
   }
   my @id = sort(@tmpid);

   #### type 1 bootstrapping will be unvalid when the sequences are quite similar
   #### say, if only 5% mutation among homologs, it will be ~5% of the mutation-involved
   #### codon columns to be selected for bootstrap, resulting the identidcal random seqs
   #### then, conversion bootstrap values will <~ 5%

   #### for type 2 bootstrapping approach, the percentages of the identical columns, mutation-
   #### involved columns and gap-involved columns in the random sequences will be equal to that
   #### in the actual sequences, which will be a prerequisition to construct random sequences.
   my (@diffpos, @idenpos, @gappos);
   my ($k1, $k2, $k3) = (0, 0, 0);
   for(my $i=0; $i<$align_len-3; $i=$i+3)
   {
      my $isgap=0;
      my $isiden=0;
      my $codon0 = $seq[0][$i].$seq[0][$i+1].$seq[0][$i+2];
      if($codon0 =~ /-/){$isgap = 1;}
      if($isgap eq 0)
      {
         for(my $j=1; $j<=$#id; $j++)
         {
            my $codon1 = $seq[$j][$i].$seq[$j][$i+1].$seq[$j][$i+2];
            if($codon1 =~ /-/){$isgap = 1; last;}
            if($codon1 eq $codon0){$isiden ++;}
         }
      }
      if($isgap eq 1){$gappos[$k1] = int($i/3); $k1++;}
      elsif($isiden eq 3){$idenpos[$k2] = int($i/3);  $k2++;}
      else{$diffpos[$k3] = int($i/3); $k3++;}
   }

   my ($count_conversion_os_Ks, $count_conversion_os_Ka) = (0, 0);
   my ($count_conversion_sb_Ks, $count_conversion_sb_Ka) = (0, 0);

   for(my $r=1; $r<=$BOOTNUM; $r++)
   {
#print "\n----------->>>".$r."\n";

       my @random_num = create_rand($align_len/3, $align_len/3);

       #my @random_num1 = create_rand2(@diffpos);
       #my @random_num2 = create_rand2(@idenpos);
       #my @random_num3 = create_rand2(@gappos);
       #my @random_num = (@random_num1, @random_num2, @random_num3);

       @random_num = permut(@random_num);

       for(my $i=0; $i<=$#random_num; $i++)
       {
          if($i eq 0)
          {
             foreach my $so ($dna_aln -> each_seq())
             {$align_seq_rand{$so -> display_id()} = "";}
          }
          for(my $j=0; $j<=$#id; $j++)
          {
             my $codon = substr($align_seq{$id[$j]}, ($random_num[$i]-1) * 3, 3);

             ##### purge off gaps after selecting the codon columns in the random sequences
             if($codon !~ /-/){$align_seq_rand{$id[$j]} .= $codon;}

#             print $random_num[$i]." ".$j." ".substr($align_seq{$id[$j]}, ($random_num[$i]-1) * 3, 3)." ";
          }
       }

       my %dna_hash;
       my $protos = Bio::SeqIO -> new(-file=> ">prot.fasta", -format=>"fasta");

       foreach my $so ($dna_aln_rand -> each_seq())
       {
          my $len = length($align_seq_rand{$so -> display_id()});
 #         print "length ...".$len."\n";
          $so -> seq($align_seq_rand{$so -> display_id()});
          $dna_hash{$so->display_id()} = $so;
          $protos -> write_seq($so -> translate());
 #         print "\n===========\n".$align_seq_rand{$so -> display_id()}."\n";
       }
       ###### without outgroup alignment
      system("clustalw -infile=prot.fasta");

      if(-z "prot.aln"){print "prot.aln is zero-sized\n"; $r=$r-1; next;}

      my $is_prot_aln = Bio::AlignIO -> new(-file=>"prot.aln", -format=>"CLUSTALW");
      my $prot_aln = $is_prot_aln -> next_aln();

      system("rm prot.fasta prot.aln prot.dnd");

      my $dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);

#      $os_check -> write_aln($dna_aln);

      my $cdsalnfile = "DNA.NUC";  #?????
      write_aln_paml($cdsalnfile, $dna_aln);

      my $seqfile = "DNA.NUC";
      my $outfile = "KV.yn00";

      set_ctl_parameters("seqfile", $seqfile,
                         "outfile", $outfile);
      system("yn00");

      my (%KaKs) = read_KYang($os_para, $sb_para, $os_sb_ortho1, $os_sb_ortho2, "KV.yn00");
      unlink($seqfile, $outfile);

      my (%Ka, %Ks);
      foreach my $pair(keys(%KaKs))
      {
         print $pair." ".$KaKs{$pair}."\n";
         my @pairvalue = split(/_/, $KaKs{$pair});
         $Ka{$pair} = $pairvalue[0];
         $Ks{$pair} = $pairvalue[1];
      }

      if($Ks{$os_para} <= $Ks{$os_sb_ortho1} && $Ks{$os_para} <= $Ks{$os_sb_ortho2})
      {
         $count_conversion_os_Ks ++;
      }

      if($Ks{$sb_para} <= $Ks{$os_sb_ortho1} && $Ks{$sb_para} <= $Ks{$os_sb_ortho2})
      {
         $count_conversion_sb_Ks ++;
      }
      
      if($Ka{$os_para} <= $Ka{$os_sb_ortho1} && $Ka{$os_para} <= $Ka{$os_sb_ortho2})
      {
         $count_conversion_os_Ka ++;
      }

      if($Ka{$sb_para} <= $Ka{$os_sb_ortho1} && $Ka{$sb_para} <= $Ka{$os_sb_ortho2})
      {
         $count_conversion_sb_Ka ++;
      }
   }
   return ($count_conversion_os_Ks/$BOOTNUM, $count_conversion_os_Ka/$BOOTNUM, 
           $count_conversion_sb_Ks/$BOOTNUM, $count_conversion_sb_Ka/$BOOTNUM);
}


sub set_ctl_parameters()
{
   my %para = @_;
   my @para = keys(%para);

   my $ctlfile = "yn00_original.ctl";
   open(CTL, $ctlfile) or die "cannot open $ctlfile due to $!\n";

   my $newctlfile = "yn00.ctl";
   open(NEW, ">".$newctlfile) or die "cannnot open $newctlfile due to $!\n";

   while(<CTL>)
   {
      my $line = $_;
      $line =~ s/^\s+//g;
      for(my $i=0; $i<=$#para; $i++)
      {
          if($line =~ /^$para[$i]/)
          {
             my $pos1 = index($_, "=");
             my $pos2 = index($_, "*");
             if($pos2 eq -1)
             {
                my $origin = substr($_, $pos1+2, length($_)-$pos1-3);
                my $update = $para{$para[$i]};
                $_ =~ s/$origin/$update/;
             }
             else
             {
                my $origin = substr($_, $pos1+2, $pos2 - $pos1);
                my $update = $para{$para[$i]}."  ";
                $_ =~ s/$origin/$update/;
             }
          }
      }
      print NEW $_;
   }

   close($newctlfile);
   close($ctlfile);
}

sub write_aln_paml()
{
        print "outputing PAML format msa ...\n";
        my($file) =$_[0];
        my($aln)  =$_[1];
        my $len  =$aln -> length();

        my @id;
        my @seq;
        my $seq_no = $aln->num_sequences;
        for(my $i=1; $i<=$seq_no; $i++)
        {
           $id[$i] = $aln -> get_seq_by_pos($i) -> display_id();
           my @tmpseq = split(//, $aln -> get_seq_by_pos($i) -> seq());
           $seq[$i] = \@tmpseq;
        }

        open(PAMLOUT, ">".$file) or die "cannot open the PAML outfile $file due to $!\n";
        print PAMLOUT "   $seq_no\t   ".$len."   I\n\n";
        for(my $i=1; $i<=$seq_no; $i++)
        { print PAMLOUT $id[$i]."\n";}
        print PAMLOUT "\n";

        my($i,$j,$startpos);
        my $gapnum = 0;

        for( $j = 0; $j <= $len/60; $j ++)
        {
                $startpos=$j*60+1;
                print PAMLOUT $startpos."\n";
                for( $i = $j*60; $i<($j+1)*60 && $i<$len; $i++ )
                {
                        print PAMLOUT $seq[1][$i];
                        if($i%3==2){print PAMLOUT " ";}
                }
                print PAMLOUT "\n";

                for(my $k = 2; $k <= $seq_no; $k ++)
                {
                   for( $i = $j*60; $i<($j+1)*60 && $i<$len; $i++ )
                   {
                        if($seq[$k][$i] ne $seq[1][$i])
                        {
                                print PAMLOUT $seq[$k][$i];
                        }
                        else
                        {
                                print PAMLOUT ".";
                        }
                        if($i%3==2){print PAMLOUT " ";}

                   }
                   print PAMLOUT "\n";
                }
        }
        close(PAMLOUT);
}

sub create_rand2()
{
   my @arr = @_;
   my $num = $#arr+1;
   my $max = $#arr+1;

   my @random_num = create_rand($num, $max);
   my @random_pos;
   for(my $i=0; $i<$num; $i++)
   {
      $random_pos[$i] = $arr[$random_num[$i]];
   }
   return @random_pos;
}

sub permut
{
   my @arr = @_;
   my $num = $#arr+1;
   my $max = $#arr+1;

   my @random_num = create_rand($num, $max);
   my @arr_permuted;
   for(my $i=0; $i<$num; $i++)
   {
      $arr_permuted[$i] = $arr[$random_num[$i]];
   }
   return @arr_permuted;
}

sub create_rand()
{
   my $num = $_[0];
   my $max = $_[1];
   my @random_num;
   for(my $i=0; $i<$num; $i++)
   {
      $random_num[$i] = int(rand($max));
   }
   return @random_num;
}

__END__