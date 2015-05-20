#!/usr/bin/perl -wi
# use strict;
use List::Util qw( min max sum );
# Cecile Ane, 2007-2011

# Requires PAUP*.

# This perl script first calls Paup to calculate the
# parsimony score of all groups of consecutive sites. This first
# step can take a very long time. Second, it calls the C++ program
# "mdl" to get the best partition (or the first k best partitions).


my $paup = "paup";   # command to call paup. Include appropriate path if needed.
my $mdl  = "~/my_mdl/src/mdl";    # command to call the C++ program mdl. Include path if needed.
my $doParsSearch = 0;
my $Ncharbase = 1000;
my $NgroupMax = 10;
my $NbestPart = 25;
my $MinBlocksInGroup = 4;
my $MaxBlocksInGroup = 12;
$excludegaps = 0;  #1=yes, 0=no.
$gap_as_character = 1;
$savetrees = 0;    #1=yes, 0=no.

$dataFile = "a.nex";
$Nsymbols = 4;     #default: DNA.

$paupFileName = "paupallgroups.nex";
$treeFile = "trees.nex";

#----- Read standard input -----------------------------------------#
foreach (@ARGV){
  if (/\bdataFile\s*=\s*([^\s]+)/i){
    $dataFile=$1;
    $root=$dataFile;
    $root =~s/.nex$//;
  }
  if (/\bdoParsSearch\s*=\s*(\w+)/i){
    if ($1 eq "yes"){ $doParsSearch=$1;}
    elsif ($1 eq "no"){ $doParsSearch=0;}
    else { die("bad option: doParsSearch = $1\n");}
  }
  if (/\bexcludegaps\s*=\s*(\w+)/i){
    if ($1 eq "yes"){   $excludegaps=1;}
    elsif ($1 eq "no"){ $excludegaps=0;}
    else { die("bad option: excludeGaps = $1\n");}
  }
  if (/\bgapaschar\s*=\s*(\w+)/i){
    if ($1 eq "yes"){   $gap_as_character=1;}
    elsif ($1 eq "no"){ $gap_as_character=0;}
    else { die("bad option: GapAsChar = $1\n");}
  }
  if (/\bncharbase\s*=\s*(\d+)/i){
    $Ncharbase=$1;
    if ($Ncharbase<1){ die("invalid NcharBase $Ncharbase");}
  }
  if (/\bngroupmax\s*=\s*(\d+)/i){ $NgroupMax=$1;}
  if (/\bnbestpart\s*=\s*(\d+)/i){ $NbestPart=$1;}
	if(/\bminblocksingroup\s*=\s*(\d+)/i){ $MinBlocksInGroup=$1; }
 if(/\bmaxblocksingroup\s*=\s*(\d+)/i){ $MaxBlocksInGroup=$1; }
}
#-------------------------------------------------------------------#

if ($excludegaps){ $gap_as_character = 0;}
#if ($excludegaps){ $extension = "_nogap";}
#elsif ($gap_as_character){ $extension = "_gapAsChar";}
#else { $extension = "_gapNotaChar";}
$extension="";
$alltreeFile = "alltrees$extension.nex";
$scoreFile = "$root.scores$extension";
$partitionScoreFile = "$root.partitions$extension";
if ($excludegaps eq 0){ $Nletters = $Nsymbols+1;}
else{                   $Nletters = $Nsymbols;}


print "dataFile is $dataFile, root is $root, doParsSearch is $doParsSearch
excludegaps is $excludegaps, gapaschar is gap_as_char, ncharbase is $Ncharbase,
ngroupmax is $NgroupMax, nbestpart is $NbestPart.\n";

#----- Read in the data file to find the Ntax and Nchar total ------#
open FHi, "<$dataFile" or die "Cannot read file $dataFile :$!";
while (<FHi>){
  if (/nchar\s*=\s*(\d+)/i){
      $Nchar = $1;
  }
  if (/ntax\s*=\s*(\d+)/i){
      $Ntax = $1;
  }
  last if ($Nchar and $Ntax);
}
close FHi;

$lastgenesize = $Nchar % $Ncharbase;
$Ngenes = ($Nchar - $lastgenesize)/$Ncharbase;
if ($lastgenesize) {$Ngenes++;}
$Ngroup = groupindex($Ngenes,$Ngenes) + 1;
print "There are $Ntax taxa, $Nchar characters, $Ngroup possible groups and\n";
print "$Ngenes consecutive small blocks of (about) equal size $Ncharbase (prior to excluding gaps).\n";

if($doParsSearch){
#---------  Make the paup file ---------------#
 unlink($scoreFile);
 unlink($treeFile);

 open FHo, ">", "$paupFileName";
 print FHo "#NEXUS\n\nset warnroot=no warntree=no warnTsave=no ";
 print FHo "increase=no maxtrees=$Ngroup monitor=no notifybeep=no;\n";
 print FHo "execute $dataFile;\n\n";
 print FHo "begin paup;\n";
 if ($gap_as_character){
     print FHo "pset gapmode=newstate;\n";
 }
 # print "NgroupMax: $NgroupMax  Ngenes: $Ngenes \n";
 foreach $startblock (1..$Ngenes){
	my $min_block = $startblock + $MinBlocksInGroup-1;
	my $max_block = min($Ngenes, $startblock + $MaxBlocksInGroup-1);
	# print "start_block, min_block, Ngenes:  [$startblock]  [$min_block] [$max_block] [$Ngenes] \n";
  foreach $endblock ($min_block..$max_block){
   # if ($endblock eq 1) last;
   # for testing purposes
   $start = startpoint($startblock);
   $end   = endpoint($endblock);

   print FHo "include $start-$end / only;\n";
   if($excludegaps){ print FHo "exclude missambig;\n";}
   print FHo "hsearch collapse=no;\n";
#	print FHo "describetrees  ;\n";
   print FHo "Pscores  1 / scorefile=$scoreFile append=yes;\n";
   if ($savetrees) {
     print FHo "savetrees from=1 to=1 file=$treeFile format=altnexus append=yes;\n";
   }
   print FHo "\n";
  }
 }

 if ($savetrees) {
  print FHo "gettrees file=$treeFile allblocks=yes;\n";
  print FHo "savetrees file=$alltreeFile format=altnexus replace=yes;\n";
 }
 print FHo "quit;\nend; [paup]";
 close FHo;

#-------- run the paup file -----------#
 system("$paup $paupFileName");

 # clean up
 # unlink($paupFileName);
} # end of: if($doParsSearch)


# print "Ngroupmax: [$NgroupMax]\n";
# print "NgroupMax+1: ", $NgroupMax+1, "\n";
my $paupout_scores = `cat $scoreFile`;
my @score_lines = split("\n", $paupout_scores);
my @score_lines_padded  =  ();
for $startblock (1..$Ngenes){
#print "NgroupMax: [$NgroupMax]\n";
	my $zzz = $startblock + $NgroupMax-1;
#print "zzz: $zzz\n";
my $min_block = $startblock + $MinBlocksInGroup-1;
my $max_block = min($Ngenes, $startblock + $MaxBlocksInGroup-1);
for $endblock ($startblock..$Ngenes){
if($endblock >= $min_block and $endblock <= $max_block){
   my $line1 = shift @score_lines;
	my $line2 = shift @score_lines;
	push @score_lines_padded, $line1;
	push @score_lines_padded, $line2;
}else{
	push @score_lines_padded, 'Tree   Length';
	push @score_lines_padded, '1    99999999';
}
}
}
my $new_score_file_name = $scoreFile . "_x";
open my $FHs, ">", "$new_score_file_name";
print $FHs (join("\n", @score_lines_padded), "\n");
close $FHs;


#-------- run the C program -----------#
# it will make the list of the best partitions of consecutive groups,
# calculate their parsimony scores, their MDLength, and sort them.
# or: it will use a heuristic to calculate the parsimony score
# and MDLength of the best partitions.

$mdl_dataoptions = "-ntax $Ntax -nchar $Nchar -scorefile $new_score_file_name -nletters $Nletters -datafile $dataFile";
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
  if ($gene<$Ngenes){ $endpoint = $gene * $Ncharbase};
  return $endpoint;
}
sub groupindex {
  my ($start,$end)=@_;
  my $ind = $start*$Ngenes - ($start*($start-1)/2) - ($Ngenes-$end) -1;
  return $ind;
}

