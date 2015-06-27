#!/usr/bin/perl -w
use strict;
use List::Util qw( min max sum );
# Cecile Ane, 2007-2011

# Requires PAUP*

# This perl script first calls Paup to calculate the
# parsimony score of all groups of consecutive sites. This first
# step can take a very long time. Second, it calls the C++ program
# "mdl" to get the best partition (or the first k best partitions).


my $paup = "paup"; # command to call paup. Include appropriate path if needed.
my $mdl  = "~/my_mdl/src/mdl"; # command to call the C++ program mdl. Include path if needed.
my $doParsSearch = 0;
my $Ncharbase = 10000;
my $NgroupMax = 1000;
my $NbestPart = 25;
my $MinBlocksInGroup = 1;
my $MaxBlocksInGroup = 12;
my $bigscore = 999999999;
my $seqFileName = 'sequences_tmp.nexus';
my $excludegaps = 0;            #1=yes, 0=no.
my $gap_as_character = 1;
my $savetrees = 0;              #1=yes, 0=no.

my $dataFile = "a.nex";
my $Nsymbols = 4;               #default: DNA.

my $paupFileName = "paupallgroups.nex";
my $treeFile = "trees.nex";

my $root;
my $extension;
my $alltreeFile;
my $scoreFile;
my $partitionScoreFile;
my $Nletters;
my $Nchar;
#my $Ntax;
my $Ntaxa;
my $lastgenesize;
my $Ngenes;
my $Ngroup;
my $startblock;
my $endblock;
my $mdl_dataoptions;
my $mdl_runoptions;


#----- Read standard input -----------------------------------------#
foreach (@ARGV) {
   if (/\bdataFile\s*=\s*([^\s]+)/i) {
      $dataFile=$1;
      $root=$dataFile;
      $root =~s/.nex$//;
   }
   if (/\bdoParsSearch\s*=\s*(\w+)/i) {
      if ($1 eq "yes") {
         $doParsSearch=$1;
      } elsif ($1 eq "no") {
         $doParsSearch=0;
      } else {
         die("bad option: doParsSearch = $1\n");
      }
   }
   if (/\bexcludegaps\s*=\s*(\w+)/i) {
      if ($1 eq "yes") {
         $excludegaps=1;
      } elsif ($1 eq "no") {
         $excludegaps=0;
      } else {
         die("bad option: excludeGaps = $1\n");
      }
   }
   if (/\bgapaschar\s*=\s*(\w+)/i) {
      if ($1 eq "yes") {
         $gap_as_character=1;
      } elsif ($1 eq "no") {
         $gap_as_character=0;
      } else {
         die("bad option: GapAsChar = $1\n");
      }
   }
   if (/\bncharbase\s*=\s*(\d+)/i) {
      $Ncharbase=$1;
      if ($Ncharbase<1) {
         die("invalid NcharBase $Ncharbase");
      }
   }
   if (/\bngroupmax\s*=\s*(\d+)/i) {
      $NgroupMax=$1;
   }
   if (/\bnbestpart\s*=\s*(\d+)/i) {
      $NbestPart=$1;
   }
   if (/\bminblocksingroup\s*=\s*(\d+)/i) {
      $MinBlocksInGroup=$1;
   }
   if (/\bmaxblocksingroup\s*=\s*(\d+)/i) {
      $MaxBlocksInGroup=$1;
   }
}
#-------------------------------------------------------------------#

if ($excludegaps) {
   $gap_as_character = 0;
}
#if ($excludegaps){ $extension = "_nogap";}
#elsif ($gap_as_character){ $extension = "_gapAsChar";}
#else { $extension = "_gapNotaChar";}
$extension="";
$alltreeFile = "alltrees$extension.nex";
$scoreFile = "$root.scores$extension";
$partitionScoreFile = "$root.partitions$extension";
if ($excludegaps eq 0) {
   $Nletters = $Nsymbols+1;
} else {
   $Nletters = $Nsymbols;
}


print "dataFile is $dataFile, root is $root, doParsSearch is $doParsSearch
excludegaps is $excludegaps, gapaschar is gap_as_char, ncharbase is $Ncharbase,
            ngroupmax is $NgroupMax, nbestpart is $NbestPart.\n";

#----- Read in the data file to find the Ntaxa and Nchar total ------#
open FHi, "<$dataFile" or die "Cannot read file $dataFile :$!";
while (<FHi>) {
   if (/nchar\s*=\s*(\d+)/i) {
      $Nchar = $1;
   }
   if (/ntax\s*=\s*(\d+)/i) {
      $Ntaxa = $1;
   }
   last if ($Nchar and $Ntaxa);
}
close FHi;

open FHi, "<$dataFile" or die "Cannot read file $dataFile :$!";
my @data_lines = <FHi>;
my @first_lines = ();
my @last_lines = ();
while (1) {
   my $line = shift @data_lines;
   push @first_lines, $line;
   last if($line =~ /^\s*matrix/);
}
my $first_lines_string = join("", @first_lines);
$first_lines_string =~ /nchar=(\d+)/;
my $Ntotalchar = $1;
# print "first lines: \n", "$first_lines_string\n";
my ($got_end, $got_semicolon) = (0, 0);
while (@data_lines) {
   my $line = pop @data_lines;
   unshift @last_lines, $line;
   $got_end = 1 if($line =~ /^\s*end/);
   $got_semicolon = 1 if($line =~ /^\s*;/);
   last if($got_end and $got_semicolon);
}

$lastgenesize = $Nchar % $Ncharbase;
$Ngenes = ($Nchar - $lastgenesize)/$Ncharbase;
if ($lastgenesize) {
   $Ngenes++;
}
$Ngroup = groupindex($Ngenes,$Ngenes) + 1;
print "There are $Ntaxa taxa, $Nchar characters, $Ngroup possible groups and\n";
print "$Ngenes consecutive small blocks of (about) equal size $Ncharbase (prior to excluding gaps).\n";

open my $FHo, ">", "$paupFileName";
print $FHo "#NEXUS\n\nset warnroot=no warntree=no warnTsave=no ";
print $FHo "increase=no maxtrees=$Ngroup monitor=no notifybeep=no;\n";
print $FHo "execute $seqFileName;\n\n";
print $FHo "begin paup;\n";
if ($gap_as_character) {
   print $FHo "pset gapmode=newstate;\n";
}
# print $FHo "include $start-$end / only;\n";
if ($excludegaps) {
   print $FHo "exclude missambig;\n";
}
print $FHo "hsearch collapse=no;\n";
#	print $FHo "describetrees  ;\n";
print $FHo "Pscores  1 / scorefile=$scoreFile append=yes;\n";
if ($savetrees) {
   print $FHo "savetrees from=1 to=1 file=$treeFile format=altnexus append=yes;\n";
}
print $FHo "\n";
if ($savetrees) {
   print $FHo "gettrees file=$treeFile allblocks=yes;\n";
   print $FHo "savetrees file=$alltreeFile format=altnexus replace=yes;\n";
}
print $FHo "quit;\nend; [paup]";
close $FHo;


my $scores_string = '';
if ($doParsSearch) {
   unlink($scoreFile);
   unlink($treeFile);
   for $startblock (1..$Ngenes) {
#  my $zzz = $startblock + $NgroupMax-1;
      my $min_block = $startblock + $MinBlocksInGroup-1;
      my $max_block = min($Ngenes, $startblock + $MaxBlocksInGroup-1);
      for $endblock ($startblock..$Ngenes) {
         if ($endblock >= $min_block and $endblock <= $max_block) { # use paup to get parsimony score & tree.
            my $first_char =  int($Ncharbase*($startblock-1));
            my $last_char = min(int($Ncharbase * $endblock - 1), $Ntotalchar-1);
            my $align_string = alignment_string($Ntaxa, \@data_lines, $first_char, $last_char);
            open my $FH1, ">", "$seqFileName" or die "Couldnt open $seqFileName for writing.\n";
            my $Nchars_in_group = $last_char - $first_char + 1;
            $first_lines_string =~ s/nchar=\d+\s*;/nchar=$Nchars_in_group;/;
            print $FH1 $first_lines_string;
            print $FH1 $align_string;
            print $FH1 join("", @last_lines);
            close $FH1;
#-------- run the paup file -----------#
            system("$paup $paupFileName");
         } else {
            $scores_string = "Tree  Length \n" . "1  $bigscore\n";
            open my $FH2, ">>", "$scoreFile" or die "couldnt open $scoreFile for writing (append).\n";
            print $FH2 $scores_string;
            close $FH2;
         }
      }
   }
# clean up
# unlink($paupFileName);
}  # end of: if($doParsSearch)


#-------- run the C program -----------#
# it will make the list of the best partitions of consecutive groups,
# calculate their parsimony scores, their MDLength, and sort them.
# or: it will use a heuristic to calculate the parsimony score
# and MDLength of the best partitions.

$mdl_dataoptions = "-ntax $Ntaxa -nchar $Nchar -scorefile $scoreFile -nletters $Nletters -datafile $dataFile";
$mdl_runoptions = "-nbestpart $NbestPart -ngroupmax $NgroupMax -o $partitionScoreFile -ncharbase $Ncharbase";
print "mdl call: $mdl $mdl_dataoptions $mdl_runoptions\n";

system("$mdl $mdl_dataoptions $mdl_runoptions");


#-------------   subroutines -----------------------------#
sub startpoint {
   my ($gene) = @_;
   my $startpoint = 1+($gene-1) * $Ncharbase;
   return $startpoint;
}
sub endpoint {
   my ($gene) = @_;
   my $endpoint = $Nchar;
   if ($gene<$Ngenes) {
      $endpoint = $gene * $Ncharbase;
   }
   ;
   return $endpoint;
}
sub groupindex {
   my ($start,$end)=@_;
   my $ind = $start*$Ngenes - ($start*($start-1)/2) - ($Ngenes-$end) -1;
   return $ind;
}


sub alignment_string{
   my $n_taxa = shift;
   my $lines = shift;       # array ref
      my $first_char = shift;  # 0-based
      my $last_char = shift;

# assume 50 chars per line (except for possibly at very end)
# also assume 2 blank lines after n_taxa lines of sequence

   my $chars_so_far = 0;
   my $first_line_index = int($first_char/50) * 7;
   my $last_line_index = (int($last_char/50) + 1) * 7 - 1;
   print "$first_char $last_char $first_line_index $last_line_index \n";
   $last_line_index = min($last_line_index, scalar @$lines - 1);
   my @selected_lines = @$lines[$first_line_index .. $last_line_index];
   return join("", @selected_lines);
}
