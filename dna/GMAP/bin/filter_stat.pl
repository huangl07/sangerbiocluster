#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$vcf);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"vcf:s"=>\$vcf,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $vcf);
open In,$fIn;
open Out,">tmp.txt";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/||/MarkerID/);
    my ($id,$type,undef)=split(/\t/,$_);
    my ($chr,$pos)=split(/\-/,$id);
    print Out join("\t",$chr,$pos),"\n";
}
close In;
close Out;
print "bcftools query -R tmp.txt -f \"%CHROM\\t%TYPE\\n\" $vcf > tmp.out";
system("bcftools query -R tmp.txt -f \"%CHROM\\t%TYPE\\n\" $vcf > tmp.out");
open In,"tmp.out";
my %stat;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($chr,$type)=split(/\t/,$_);
    if($type =~ /INDEL/){$type = "INDEL"}
    $stat{$chr}{$type}++;
}
close In;
open Out,">$fOut";
print Out "CHROM\tSNP\tINDEL\n";
foreach my $chr(sort keys %stat){
    $stat{$chr}{SNP}||=0;
    $stat{$chr}{INDEL}||=0;
    print Out $chr,"\t",$stat{$chr}{SNP},"\t",$stat{$chr}{INDEL},"\n"
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
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"vcf:s"=>\$vcf,
  -h         Help

USAGE
        print $usage;
        exit;
}
