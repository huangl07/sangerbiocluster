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
my ($fIn,$fOut,$groupfile,$wild,$mutate);
GetOptions(
                "help|?" =>\&USAGE,
                "vcf:s"=>\$fIn,
                "out:s"=>\$fOut,
                "group:s"=>\$groupfile
                                ) or &USAGE;
&USAGE unless ($fIn and $fOut);

open In,$groupfile;
my ($sampleID,$group,$low,$high);
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	($sampleID,$group,$low,$high)=split(/\s+/,$_) ;
	if($group eq "P1"){
		$wild=$sampleID;
	}elsif($group eq "B2"){
		$mutate=$sampleID;
	}
}
$wild||="-";
$mutate||="-";
close In;

open In,"zcat $fIn|";
open Out,">$fOut";
my @Indi;
my %samples;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/ || /^##/);
    if(/^#/){
        my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples)=split(/\t/,$_);
        for(my $i=0;$i<@samples;$i++){
            push @Indi,$samples[$i];
            $samples{$samples[$i]}=1;
        }
        if(!exists $samples{$wild}){die "error in wild parent: $wild\n"}
        if(!exists $samples{$mutate}){die "error in wild parent: $mutate\n"}
        print Out join("\t","#CHROM","POS","POS","REF","ALT",$wild,$mutate),"\n";
    }else{
        my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples)=split(/\t/,$_);
        my @format=split(/\:/,$format);
        my $wildtype;
        my $mutatetype;
        my @GT=split(/,/,join(",",$ref,$alt));
        for(my $i=0;$i<@samples;$i++){
            next if($Indi[$i] ne $wild && $Indi[$i] ne $mutate);
            my @info=split(/\:/,$samples[$i]);
            for(my $j=0;$j<@info;$j++){
                if($format[$j] eq "GT"){
                    if($info[$j] ne "./."){
                        my ($g1,$g2)=split("/",$info[$j]);
                        $wildtype=join("",sort($GT[$g1],$GT[$g2])) if($Indi[$i] eq $wild);
                        $mutatetype=join("",sort($GT[$g1],$GT[$g2])) if($Indi[$i] eq $mutate);
                        last;
                    }else{
                        $wildtype="NN" if($Indi[$i] eq $wild);
                        $mutatetype="NN" if($Indi[$i] eq $mutate);
                    }
                }
            }
        }
        if(($wildtype eq "GG" && $mutatetype eq "AA") ||($wildtype eq "GG" && $mutatetype eq "AG")|| ($wildtype eq "CC" && $mutatetype eq "TT") ||($wildtype eq "CC" && $mutatetype eq "CT")){
            print Out join("\t",$chrom,$pos,$pos,$ref,$alt,$wildtype,$mutatetype),"\n";
        }
    }

}
close Out;
close In;
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
    -vcf        vcf file
    -group      group file
    -out        output file

USAGE
        print $usage;
        exit;
}