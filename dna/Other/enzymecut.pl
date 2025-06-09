#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fa,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
	"output:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fa and $fOut);
my $MseI="AATT";
my $TaqaI="TCGA";
my @enzyme;
push @enzyme,$MseI;
push @enzyme,$TaqaI;
open In,$fa;
open Out,">$fOut";
print Out join("\t","#chr","pos1","pos2","length","enzyme1","enzyme2"),"\n";
$/=">";
my %stat;
my %estat;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($id,@line)=split(/\n/,$_);
    $id=(split(/\s+/,$id))[0];
    my $seq=join("",@line);
    my $rseq=reverse($seq);
    $rseq=~tr/ATGC/TACG/;
    my %pos;
    foreach my $enzyme(@enzyme){
        while($seq=~/$enzyme/ig){
            my $pos=pos($seq);
            $pos{$pos}=$enzyme;

        }
        while($rseq=~/$enzyme/ig){
            my $pos=pos($rseq);
            $pos{$pos}=$enzyme;
        }
    }
    my @pos=sort {$a <=> $b} keys %pos;
    for(my $i=0;$i<@pos-1;$i++){
        if($pos{$pos[$i]} ne $pos{$pos[$i+1]}){
            my $lenth=$pos[$i+1]-$pos[$i];
            print Out join("\t",$id,$pos[$i],$pos[$i+1],$lenth,$pos{$pos[$i]},$pos{$pos[$i+1]}),"\n";
            $stat{"250~300"}++ if( $lenth >= 250 && $lenth < 300);
            $stat{"300~350"}++ if( $lenth >= 300 && $lenth < 350);
            $stat{"350~400"}++ if( $lenth >= 350 && $lenth < 400);
            $stat{"400~450"}++ if( $lenth >= 400 && $lenth < 450);
            $stat{"450~500"}++ if( $lenth >= 450 && $lenth < 500);
            $stat{"500~550"}++ if( $lenth >= 500 && $lenth < 550);
            $stat{"550~600"}++ if( $lenth >= 550 && $lenth < 600);

            $estat{"250~300"}++ if( $lenth >= 250-120 && $lenth < 300-120);
            $estat{"300~350"}++ if( $lenth >= 300-120 && $lenth < 350-120);
            $estat{"350~400"}++ if( $lenth >= 350-120 && $lenth < 400-120);
            $estat{"400~450"}++ if( $lenth >= 400-120 && $lenth < 450-120);
            $estat{"450~500"}++ if( $lenth >= 450-120 && $lenth < 500-120);
            $estat{"500~550"}++ if( $lenth >= 500-120 && $lenth < 550-120);
            $estat{"550~600"}++ if( $lenth >= 550-120 && $lenth < 600-120);
        }
    }
}
close In;
close Out;
open Out,">$fOut.withoutPrimer.stat";
foreach my $key (sort keys %stat){
    print Out join("\t",$key,$stat{$key}),"\n";
}
close Out;
open Out,">$fOut.stat";
foreach my $key (sort keys %estat){
    print Out join("\t",$key,$estat{$key}),"\n";
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
	perl $Script -fa -output

Usage:
  Options:

    "fa:s"=>\$fa,
	"output:s"=>\$fOut,

USAGE
        print $usage;
        exit;
}
