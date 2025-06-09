#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamstat,$metric,$fOut,$insert,$depth,$Key,$fai);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"b:s"=>\$bamstat,
	"i:s"=>\$insert,
	"c:s"=>\$depth,
	"o:s"=>\$fOut,
	"k:s"=>\$Key,
	"f:s"=>\$fai
			) or &USAGE;
&USAGE unless ($fOut and $bamstat and $insert and $depth);
open In,$bamstat;
open IS,">$insert";
open COV,">$depth";
my %mapstat;
my %cov;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/ );
	if (/^SN/) {
		if (/raw total sequences:/) {
			$mapstat{total}=(split(/\t/,$_))[2];
		}
		if (/reads mapped\:/) {
			$mapstat{mapped}=(split(/\t/,$_))[2];
		}
		if (/reads mapped and paired:/) {
			$mapstat{paired}=(split(/\t/,$_))[2];
		}
		if (/reads properly paired:/) {
			$mapstat{properly}=(split(/\t/,$_))[2];
		}
		if (/reads duplicated:/) {
			$mapstat{dup}=(split(/\t/,$_))[2];
		}
		if (/insert size average:/) {
			$mapstat{insert}=(split(/\t/,$_))[2];
		}
		if (/bases mapped:/) {
			$mapstat{coverage}=(split(/\t/,$_))[2];
		}
	}
	if (/^IS/) {
		my (undef,$insert,$depth,undef,undef)=split("\t",$_);
		print IS $insert,"\t",$depth,"\n";
	}
	if (/^COV/) {
		my (undef,undef,$deps,$cova)=split("\t",$_);
		$cov{1}+=$cova if($deps >= 1);
		$cov{4}+=$cova if($deps >= 4);
		print COV $deps,"\t",$cova,"\n";
	}
}
close IS;
close COV;
close In;

open FAI,$fai;
my $totallen=0;
while (<FAI>) {
	chomp;
	my (undef,$len,undef,undef,undef)=split("\t",$_);
	$totallen+=$len;
}
close FAI;

open Out,">$fOut";
print Out "Sample ID\t","Mapped Ratio(%)\t","Proper Ratio(%)\t","Duplicate Ratio(%)\t","Average Insert Size\t","Average Depth\t","Real Depth\t","Coverage(%) (>=1x)\t","Coverage(%) (>=4x)\n";
print Out $Key,"\t";
print Out sprintf("%.2f",$mapstat{mapped}/$mapstat{total}*100),"\t";
print Out sprintf("%.2f",$mapstat{properly}/$mapstat{total}*100),"\t";
print Out sprintf("%.2f",$mapstat{dup}/$mapstat{total}*100),"\t";
print Out $mapstat{insert},"\t";
print Out sprintf("%.2f",$mapstat{coverage}/$totallen),"\t";
print Out sprintf("%.2f",$mapstat{coverage}/$cov{1}),"\t";
print Out sprintf("%.2f",$cov{1}/$totallen*100),"\t";
print Out sprintf("%.2f",$cov{4}/$totallen*100),"\n";
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        yiwei.tang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -b a.mapstat -i a.insert -c a.depth -o a.result.stat -k sample_id

Usage:
  Options:
  -b	<file>	samtools stat file
  -i	<file>	output insert file
  -c	<file>	output cover depth file
  -o	<file>	output stat file
  -k	<str>	sample name
  -f	<file>	reference sequences index
  -h         Help

USAGE
        print $usage;
        exit;
}
