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
my ($fIn,$fOut,$group,$P1,$P2,$B1,$B2,$popt,$minindex);
GetOptions(
				"help|?" =>\&USAGE,
				"table:s"=>\$fIn,
				"out:s"=>\$fOut,
				"group:s"=>\$group,
				"popt:s"=>\$popt,
				"minindex:f"=>\$minindex
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$group;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sampleID,$group,$low,$high)=split(/\s+/,$_) ;
	if($group eq "P1"){
		$P1=$sampleID;
	}elsif($group eq "P2"){
		$P2=$sampleID;
	}elsif($group eq "B1"){
		$B1=$sampleID;
	}elsif($group eq "B2"){
		$B2=$sampleID;
	}
}
$B1||="-";
$P1||="-";
$P2||="-";
$B2||="-";
close In;
$popt ||="F2";
open In,$fIn;
open Out,">$fOut";
if ($B1 eq "-") {
	print Out "#chr\tpos\tdelta\tdepth\n";
}else{
	print Out "#chr\tpos\tindex1\tindex2\tdelta\tn1\tn2\tn3\tn4\n";
}
my @head;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/ ||/^CHROM/) {
		my ($chr,$pos,$ref,$alt,$type,@sample)=split(/\s+/,$_);
		push @head,@sample;
		pop(@head);
	}else{
		my ($chr,$pos,$ref,$alt,$type,@sample)=split(/\s+/,$_);
		my %info;
		for (my $i=0;$i<@head;$i++) {
			my ($id,$info)=(split(/\./,$head[$i]));
			$info{$id}{$info}=$sample[$i];
		}
		my %ad;
		foreach my $sample (keys %info) {
			my @ad=split(/\,/,$info{$sample}{AD});
			for (my $i=0;$i<@ad;$i++) {
				$ad{$sample}{$ref}=$ad[0];
				$ad{$sample}{$alt}=$ad[1];
			}
		}
		if ($popt eq "F2") {
			if ($B1 eq "-") {
				my $index;
				if ($P1 ne "-" && $P2 eq "-") {
					my @gtP=split(/\//,$info{$P1}{GT});
					if ($gtP[0] eq $ref) {
						$index=$ad{$B2}{$alt}/$info{$B2}{DP};
					}else{
						$index=$ad{$B2}{$ref}/$info{$B2}{DP};
					}
				}elsif ($P2 ne "-" && $P1 eq "-") {
					my @gtP=split(/\//,$info{$P2}{GT});
					$index=$ad{$B2}{$gtP[0]}/$info{$B2}{DP};
				}elsif ($P1 ne "-" && $P2 ne "-") {
					my @gtP=split(/\//,$info{$P2}{GT});
					$index=$ad{$B2}{$gtP[0]}/$info{$B2}{DP};
				}elsif ($P2 eq "-" && $P1 eq "-") {
					$index=$ad{$B2}{$alt}/$info{$B2}{DP};
				}
				next if($index<$minindex);
				print Out join("\t",$chr,$pos,$index,$info{$B2}{DP}),"\n";
			}else{
				my ($index1,$index2,$delta);
				if ($P1 ne "-" && $P2 eq "-") {
					my @gtP=split(/\//,$info{$P1}{GT});
					my $numerator=$ref;
					$numerator=$alt if($gtP[0] eq $ref);
					if(!exists $ad{$B2}{$numerator}){
						print $_,"\n";
						print $numerator,"\n";
						print Dumper %ad;
						print Dumper %info;
						die;
					};
					$index1=$ad{$B2}{$numerator}/$info{$B2}{DP};
					$index2=$ad{$B1}{$numerator}/$info{$B1}{DP};
					$delta=$index1-$index2;
				}elsif ($P2 ne "-" && $P1 eq "-") {
					my @gtP=split(/\//,$info{$P2}{GT});
					my $numerator=$gtP[0];
					$index1=$ad{$B2}{$numerator}/$info{$B2}{DP};
					$index2=$ad{$B1}{$numerator}/$info{$B1}{DP};
					$delta=$index1-$index2;
				}elsif ($P1 ne "-" && $P2 ne "-") {
					my @gtP=split(/\//,$info{$P2}{GT});
					my $numerator=$gtP[0];
					$index1=$ad{$B2}{$numerator}/$info{$B2}{DP};
					$index2=$ad{$B1}{$numerator}/$info{$B1}{DP};
					$delta=$index1-$index2;
				}elsif ($P1 eq "-" && $P2 eq "-") {
					$index1=$ad{$B2}{$alt}/$info{$B2}{DP};
					$index2=$ad{$B1}{$alt}/$info{$B1}{DP};
					$delta=abs($index1-$index2);
				}
				print Out join("\t",$chr,$pos,$index1,$index2,$delta,$ad{$B1}{$ref},$ad{$B2}{$ref},$ad{$B1}{$alt},$ad{$B2}{$alt}),"\n";
			}
		}else{#F1
			if ($B1 eq "-") {
				my $index=$ad{$B2}{$alt}/$info{$B2}{DP};
				print Out join("\t",$chr,$pos,$index,$info{$B2}{DP}),"\n";
			}else{
				my ($index1,$index2,$delta);
				$index1=$ad{$B2}{$alt}/$info{$B2}{DP};
				$index2=$ad{$B1}{$alt}/$info{$B1}{DP};
				$delta=abs($index1-$index2);
				print Out join("\t",$chr,$pos,$index1,$index2,$delta,$ad{$B1}{$ref},$ad{$B2}{$ref},$ad{$B1}{$alt},$ad{$B2}{$alt}),"\n";
			}
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

	-table input table for calculate
	-out	output filename for result
	--group	<file>	input group file
		mbid must be given

		wpid    mpid    phdep   pldep
        wbid    mbip    bhdep   bldep

	-popt	population type

USAGE
	print $usage;
	exit;
}
