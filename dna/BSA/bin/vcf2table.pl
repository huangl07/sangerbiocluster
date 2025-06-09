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
my ($fIn,$fOut,$P1,$P2,$B1,$B2,$Pldep,$Bldep,$Bhdep,$Phdep,$popt,$Vtype,$group);
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$fIn,
	"out:s"=>\$fOut,
	"group:s"=>\$group,
	"popt:s"=>\$popt,
	"vtype:s"=>\$Vtype
				) or &USAGE;
&USAGE unless ($fIn and $fOut );
$Pldep||=10;
$Phdep||=1000;
$Bldep||=10;
$Bhdep||=1000;
$popt||="F2";
$Vtype||="ALL";
$P1="-";
$P2="-";
$B1="-";
$B2="-";
open In,$group;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sampleID,$group,$high,$low)=split(/\s+/,$_);
	if($group eq "P1"){
		$P1=$sampleID;
		$Pldep=$low;
		$Phdep=$high;
	}elsif($group eq "P2"){
		$P2=$sampleID;
		$Pldep=$low;
		$Phdep=$high;
	}elsif($group eq "B1"){
		$B1=$sampleID;
		$Bldep=$low;
		$Bhdep=$high;
	}elsif($group eq "B2"){
		$B2=$sampleID;
		$Bldep=$low;
		$Bhdep=$high;
	}
}
close In;
if(!defined $B1 && !defined $B2){
	die "error group file die $B1 $B2";
}elsif(!defined $B2){
	die "$B2";
}


if ($fIn =~ /\.gz$/){
	open In,"gzip -dc $fIn|";
}else{
	open In,$fIn;
}
open Out,">$fOut";
my $bed=$fOut;
$bed=~s/table/anno/g;
open BED,">$bed";
my @Sample;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\t/,$_);
		push @Sample,@info;
		my %indi;
		foreach my $s(@info){
			$indi{$s}=1;
		}
		my @out;
		die $P1 if($P1 ne "-" && !exists $indi{$P1});
		die $P2 if($P2 ne "-" && !exists $indi{$P2});
		die $B1 if($B1 ne "-" && !exists $indi{$B1});
		die $B2 if($B2 ne "-" && !exists $indi{$B2});
		push @out,join("\t","$P1.GT","$P1.AD","$P1.DP") if (defined $P1 && $P1 ne "-");
		push @out,join("\t","$P2.GT","$P2.AD","$P2.DP") if (defined $P2 && $P2 ne "-");
		push @out,join("\t","$B1.GT","$B1.AD","$B1.DP") if (defined $B1 && $B1 ne "-");
		push @out,join("\t","$B2.GT","$B2.AD","$B2.DP") if (defined $B2 && $B2 ne "-");
		print BED join("\t","CHROM\tPOS\tRef\tAlt\tVtype",@out,"Allele\tAnnotation\tPutativeImpact\tGeneName\tGeneID\tFeatureType\tFeatureID\tTranscriptBioType\tRank\tHGVS.c\tHVGS.p\tcDNA_position_len\tCDS_position_len\tProtein_position_len\tDistance\tERRORS"),"\n";
		print Out join("\t","CHROM\tPOS\tRef\tAlt\tVtype",@out,"ANN"),"\n";
	}else{
		my (%samples,$vtype,@ann,$baseinfo,$allenum);
		$vtype=readvcf($_,\%samples,\@Sample,\@ann,\$baseinfo,\$allenum);
		if ($vtype ne $Vtype && $Vtype ne "ALL"){
			$stat{vtype}++;
			next;
		}
		if (defined $B1 && $B1 ne "-" && $samples{$B1}{GT} ne "NN" && ($samples{$B1}{DP} <= $Bldep|| $samples{$B1}{DP} >= $Bhdep)){
			$stat{B1depth}++;
			next;
		}
		if ($samples{$B2}{GT} ne "NN" && ($samples{$B2}{DP} <= $Bldep || $samples{$B2}{DP} >= $Bhdep)){
			$stat{B2depth}++;
			next;
		}
		if ($P1 ne "-" && $samples{$P1}{GT} ne "NN"&&($samples{$P1}{DP} <= $Pldep || $samples{$P1}{DP} >= $Phdep)){
			$stat{P1depth}++;
			next;
		}
		if ($P2 ne "-" && $samples{$P2}{GT} ne "NN"&&($samples{$P2}{DP} <= $Pldep || $samples{$P2}{DP} >= $Phdep)){
			$stat{P2depth}++;
			next;
		}
		if ($B1 ne "-" && ($samples{$B1}{GT} eq "NN" )){
			$stat{B1missing}++;
			next;
		}
		if ($samples{$B2}{GT} eq "NN" ){
			$stat{B2missing}++;
			next;
		}
		if ($P1 ne "-" && ($samples{$P1}{GT} eq "NN")){
			$stat{P1missing}++;
			next;
		}
		if ($P2 ne "-" && ($samples{$P2}{GT} eq "NN")){
			$stat{P2missing}++;
			next;
		}
		if ($allenum > 2){
			$stat{allenumMoreThan2}++;
			next;
		};
		if ($popt ne "F1" && $P1 ne "-" && $P2 ne "-" &&defined $P1 && defined $P2 && $samples{$P1}{GT} eq $samples{$P2}{GT}){
			$stat{PMsameGenotype}++;
			next;
		};
		my @g1=split(/\||\//,$samples{$P1}{GT}) if($P1 ne "-" && defined $P1);
		my @g2=split(/\||\//,$samples{$P2}{GT}) if ($P2 ne "-" && defined $P2);
		if ($popt ne "F1" && ((defined $P1 && $P1 ne "-" && $g1[0] ne $g1[1]) || (defined $P2 && $P2 ne "-" && $g2[0] ne $g2[1]))){
			$stat{popErrF2buthete}++;
			next;
		};
		if ($popt eq "F1" && ((defined $P1 && $P1 ne "-" && $g1[0] eq $g1[1]) ||(defined $P2 && $P2 ne "-"  && $g2[0] eq $g2[1]))){
			$stat{popErrF1buthomo}++;
			next;
		}
		my ($ge1,$ge2)=split("/",$samples{$B2}{GT});
		if($popt eq "F2" && $ge1 eq $ge2 && $B1 ne "-" && $samples{$B2}{GT} eq $samples{$B1}{GT}){
			$stat{popF2BulkGTsame}++;
			next;
		}
		if($popt eq "F2" && $P1 ne "-" && $samples{$B2}{GT} eq $samples{$P1}{GT}){
			$stat{popF2B2sameP1}++;
			next;
		}
		$stat{Pass}++;
		my @out;
		push @out,join("\t",$samples{$P1}{GT},$samples{$P1}{AD},$samples{$P1}{DP}) if (defined $P1 && $P1 ne "-");
		push @out,join("\t",$samples{$P2}{GT},$samples{$P2}{AD},$samples{$P2}{DP}) if (defined $P2 && $P2 ne "-");
		push @out,join("\t",$samples{$B1}{GT},$samples{$B1}{AD},$samples{$B1}{DP}) if (defined $B1 && $B1 ne "-");
		push @out,join("\t",$samples{$B2}{GT},$samples{$B2}{AD},$samples{$B2}{DP}) if (defined $B2 && $B2 ne "-");
		print Out join("\t",$baseinfo,$vtype,@out,join("",@ann)),"\n";
		foreach my $a(@ann){
			$a=~s/\|/\t/g;
			print BED join("\t",$baseinfo,$vtype,@out,$a),"\n";
		}
	}
}
close In;
close Out;
close BED;
open Out,">$fOut.stat";
foreach my $key (sort keys %stat){
	print Out $key,"\t",$stat{$key},"\n";
}
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub readvcf{
	my ($line,$sample,$Sample,$anninfo,$baseinfo,$allenum)=@_;
	my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\s+/,$line);
	$$baseinfo=join("\t",$chrom,$pos,$ref,$alt);
	my @alles=split(",",join(",",$ref,$alt));
	$$allenum=scalar @alles;
	if($info=~/ANN=([^\;]*)/g){
		my @ann=split(/\,/,$1);
		for (my $i=0;$i<@ann;$i++) {
			my @str=split(/\|/,$ann[$i]);
			for (my $j=0;$j<16;$j++){
				$str[$j]="--" if(!defined $str[$j]|| $str[$j] eq "");
			}
			my $ann=join("|",@str);
			push @{$anninfo},$ann;
		}
	}
	my %len;
	for (my $i=0;$i<@alles;$i++) {
		$len{length($alles[$i])}=1;
	}
	my $type="SNP";
	if (scalar keys %len > 1) {
		$type="INDEL";
	}
	my @format=split(/\:/,$format);
	for (my $i=0;$i<@info;$i++) {
		my @infos=split(/\:/,$info[$i]);
		for (my $j=0;$j<@infos;$j++) {
			if ($format[$j] eq "GT") {
				if ($infos[$j] =~ /\./) {
					$$sample{$$Sample[$i]}{$format[$j]}="NN";
					$$sample{$$Sample[$i]}{DP}=0;
					$$sample{$$Sample[$i]}{AD}=0;
				}else{
					my @gt=split(/\/|\//,$infos[$j]);
					$$sample{$$Sample[$i]}{$format[$j]}=join("/",sort($alles[$gt[0]],$alles[$gt[1]]));
				}
			}
			if ($format[$j] eq "AD") {
				$$sample{$$Sample[$i]}{$format[$j]}=$infos[$j];
			}
			if ($format[$j] eq "DP") {
				$$sample{$$Sample[$i]}{$format[$j]}=$infos[$j];
			}
		}
	}
	return $type;
}
sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-vcf	<file>	input file 
	-out	<file>	output file
	--group	<file>	input group file
		mbid must be given

		wpid    P1	phdep   pldep
		wbid    P2	mbip    bhdep  

	-popt	<str>	population type default F2
	-vtype	<str>	varaint type default ALL
USAGE
	print $usage;
	exit;
}
