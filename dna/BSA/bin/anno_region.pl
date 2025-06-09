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
my ($fIn,$fOut,$Region,$P1,$P2,$B1,$B2,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"table:s"=>\$fIn,
				"out:s"=>\$fOut,
				"region:s"=>\$Region,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$Region;
my %region;
while(<In>){
	chomp;
	next if(/^$/ || $_ eq "");
	my ($chr,$start,$end,undef)=split;
	$region{$chr}{$start}=$end;
}
close In;
open In,$fIn;
open Out,">$fOut";
while(<In>){
	chomp;
	next if(/^$/ || $_ eq "");
	if(/^CHROM/){
		print Out "REGION\t$_\n";
		next;
	}
	my ($chr,$pos,$info)=split(/\s+/,$_,3);
	foreach my $start(sort keys %{$region{$chr}}){
		if($pos >= $start && $pos <= $region{$chr}{$start}){
			print Out join(":",$chr,$start,$region{$chr}{$start}),"\t$_\n";
		}
	}
}
close In;
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"111s\n";
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
