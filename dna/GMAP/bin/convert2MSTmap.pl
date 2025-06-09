#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$popt,$homo,$seg,$mis,$chr,$marker);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"popt:s"=>\$popt,
	"chr:s"=>\$chr,
    "marker:s"=>\$marker,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $popt and $chr);

open In,$marker;
my @indi;
while(<In>){
    chomp;
    next if(/^$/||$_ eq "");
    (undef,undef,@indi)=split(/\s+/,$_);
    last;
}
close In;
open In,$fIn;
my @start;
my @end;
my %genotype;
my @center;
while (<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    if(/bin start/){
        (undef,@start)=split(/,/,$_);
    }elsif(/bin end/){
        (undef,@end)=split(/,/,$_);
    }elsif(/bin center/){
        (undef,@center)=split(/,/,$_);
    }else{
        my ($id,@info)=split(/,/,$_);
        for(my $i=0;$i<@info;$i++){
            if($info[$i] eq "h" ){
                $info[$i] = "X";
            }
            if($info[$i] eq "X" && $popt eq "DH"){
                $info[$i] = "-";
            }
            $genotype{$center[$i]}{id}=$start[$i]."_".$end[$i];
            $genotype{$center[$i]}{type}{$id}=$info[$i];
        }
    }
}
close In;

my $nloc=scalar keys %genotype;
my $id="LG";
my $nind=scalar @indi;
open Out,">$fOut";
print Out join("\t","MarkerID",@indi),"\n";
foreach my $keys (sort {$a<=>$b} keys %genotype){
    my @out;
    foreach my $indi(@indi){
        $genotype{$keys}{type}{$indi}||="U";
        push @out,$genotype{$keys}{type}{$indi};
    }
    print Out join("\t","$chr\-".$genotype{$keys}{id},@out),"\n";
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
  -popt     <str>	population type "RIL\\d | DH"
  -chr      <chr>   input keys of chr for bins
  -h         Help

USAGE
        print $usage;
        exit;
}
