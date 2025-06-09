#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$SSR);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"SSR:s"=>\$SSR,
    "output:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $SSR);
open In,$fIn;
open Out,">$fOut";
my @indi;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    if(/MarkerID/){
        (undef,undef,@indi)=split(/\s+/,$_);
        print Out $_,"\n";
    }else{
        print Out $_,"\n";
    }
}
close In;
open In,$SSR;
my @indis;
while(<In>){
    chomp;
    next if($_ eq ""||/^$/);
    if(/MarkerID/){
        (undef,undef,@indis)=split(/\t/,$_);
        next;
    }
    my ($markerID,$type,@info)=split(/\t/,$_);
    my %genotype;
    for(my $i=0;$i<@info;$i++){
        $info[$i]="--" if($info[$i] eq "");
        $genotype{$indis[$i]}=$info[$i];
    }
    my @out;
    push @out,$markerID;
    push @out,$type;
    for(my $i=0;$i<@indi;$i++){
        if(!exists $genotype{$indi[$i]}){
            print "\n\n";
            print Dumper keys %genotype;
            print $indi[$i];die;
        }
        push @out,$genotype{$indi[$i]};
    }
    print Out join("\t",@out),"\n";
}
close In;
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
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"SSR:s"=>\$vcf,
  -h         Help

USAGE
        print $usage;
        exit;
}
