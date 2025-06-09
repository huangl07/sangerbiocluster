#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$fanno,$P1,$P2,$B1,$B2,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"table:s"=>\$fIn,
				"out:s"=>\$fOut,
				"anno:s"=>\$fanno,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fanno;
my %info;
my $ninfo;
my @ID;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/ );
    if(/GeneID/){
        (undef,@ID)=split(/\t/,$_);
        $ninfo=scalar @ID;
        next;
    }
    my ($pos,@info)=split(/\t/,$_);
    for (my $i=0;$i<$ninfo;$i++){
        $info[$i]="--" if(!defined $info[$i] || $info[$i] eq "");
    }
    my @mtype=split(/:/,$pos);
    $info{$mtype[1]}=join("\t",@info);
    
}
close In;
open In,$fIn;
open Out,">$fOut";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    if(/CHROM/){
        print Out $_,"\t",join("\t",@ID),"\n";
        next;
    }
    my @a=split(/\t/,$_);
    if(exists $info{$a[-10]}){
        print Out join("\t",@a,$info{$a[-10]}),"\n";
    }else{
        my @out;
        for(my $i=0;$i<$ninfo;$i++){
            push @out,"--";
        }
        print Out join("\t",@a,@out),"\n";
    }
}
close In;
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:

				"table:s"=>\$fIn,
				"out:s"=>\$fOut,
				"anno:s"=>\$anno,
	-popt	population type

USAGE
	print $usage;
	exit;
}
