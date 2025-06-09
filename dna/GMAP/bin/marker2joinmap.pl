#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin qw ($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($fIn,$fOut,$map);
GetOptions(
    "help|?"=>\&USAGE,
    "input:s"=>\$fIn,
    "output:s"=>\$fOut,
    "map:s"=>\$map
) or &USAGE;
&USAGE unless ($fIn and $fOut);
open IN,$fIn;
my $nindi=0;
my $nloc=0;
my @marker;
while(<IN>){
    chomp;
    next if ($_ eq "" || /^$/ ||/^#/);
    my ($markerid,$type,$phase,@geno)=split(/\s+/,$_);
    push @marker,join("\t",$markerid,$type,$phase,@geno);    
    $nloc++;
    $nindi=scalar @geno;
}
close IN;
open Out,">$fOut.loc";
print Out "nind=$nindi\n";
print Out "nloc=$nloc\n";
print Out "popt=CP\n";
print Out "name=name\n";
print Out join("\n",@marker);
close Out;
open Out,">$fOut.trt";
print Out "nind=$nindi\n";
print Out "ntrt=1\n";
print Out "miss=*\n";
print Out "trt\n";
for(my $i=0;$i<$nindi;$i++){
    print Out $i,"\n";
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -input	<file>	input file name
  -output	<file>	input keys of file name

  -h         Help

USAGE
        print $usage;
        exit;
}
