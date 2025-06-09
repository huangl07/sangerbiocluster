#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($rawdup , $cleandup , $qc , $dOut,$pollution,$add);
GetOptions(
				"help|?" =>\&USAGE,
				"pollution:s"=>\$pollution,
				"qc:s"=>\$qc,
				"o:s"=>\$dOut,
				"add"=>\$add,
				) or &USAGE;
&USAGE unless ($qc and $dOut);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
my %stat;
open In,"$qc/stat.list";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$stats)=split(/\s+/,$_);
	open Stat,$stats;
	while (<Stat>) {
		next if ($_ eq "" || /^$/|| /#/);
		($sample,$stat{$sample}{rawread},$stat{$sample}{rawbase},$stat{$sample}{rawq20},$stat{$sample}{rawq30},$stat{$sample}{rawGC},$stat{$sample}{ada},$stat{$sample}{cleanread},$stat{$sample}{cleanbase},$stat{$sample}{cleanq20},$stat{$sample}{cleanq30},$stat{$sample}{cleanGC},$stat{$sample}{dup})=split(/\t/,$_);
#		next if ($_ eq "" || /^$/|| /#/);##sampleID	rawreads	rawdata	rawq20	rawq30	rawGC	A	T	C	G	N	LowQ reads	N reads	ada	cleanreads	cleandata	cleanq20	cleanq30	cleanGC	A	T	C	G	N	dup
#		($sample,$stat{$sample}{rawread},$stat{$sample}{rawbase},$stat{$sample}{rawq20},$stat{$sample}{rawq30},$stat{$sample}{rawGC},$stat{$sample}{rA},$stat{$sample}{rT},$stat{$sample}{rC},$stat{$sample}{rG},$stat{$sample}{rN},$stat{$sample}{lowreads},$stat{$sample}{Nreads},$stat{$sample}{ada},$stat{$sample}{cleanread},$stat{$sample}{cleanbase},$stat{$sample}{cleanq20},$stat{$sample}{cleanq30},$stat{$sample}{cleanGC},$stat{$sample}{cA},$stat{$sample}{cT},$stat{$sample}{cC},$stat{$sample}{cG},$stat{$sample}{cN},$stat{$sample}{dup})=split(/\t/,$_);
	}
	close Stat;
}
close In;
#open In,$pollution;
#while (<In>) {
#	chomp;
#	next if ($_ eq "" ||/^$/);
#	my ($sample,$po1,$po2,$po3)=split(/\s+/,$_);
#	$stat{$sample}{pollutions}=$po1;
#	$stat{$sample}{pollutionr}=$po2;
#	$stat{$sample}{pollutionn}=$po3;
#}
close In;
open Out,">$dOut/qc-report.xls";

print Out join("\t","#SampleID","Raw Reads","Raw Base","Raw Q20%","Raw Q30%","Raw GC%","Clean Reads","Clean Base","Clean Q20%","Clean Q30","Clean GC%","Dup%","Adapter%"),"\n";#,"Pollution-species","Pollution-rate%","Pollution-number"),"\n";
foreach my $sample (sort keys %stat) {
	#my $samples=$sample;
	print Out join("\t",$sample,$stat{$sample}{rawread}/2,$stat{$sample}{rawbase},sprintf("%.2f",100*$stat{$sample}{rawq20}),sprintf("%.2f",100*$stat{$sample}{rawq30}),sprintf("%.2f",100*$stat{$sample}{rawGC}),$stat{$sample}{cleanread}/2,$stat{$sample}{cleanbase},sprintf("%.2f",100*$stat{$sample}{cleanq20}),sprintf("%.2f",100*$stat{$sample}{cleanq30}),sprintf("%.2f",100*$stat{$sample}{cleanGC}),sprintf("%.2f",100*$stat{$sample}{dup})),"\n";#,sprintf("%.2f",100*$stat{$sample}{ada}),$stat{$sample}{pollutions},$stat{$sample}{pollutionr},$stat{$sample}{pollutionn}),"\n";
}
close Out;
mkdir "$dOut/01_RawdataStat" if (!-d "$dOut/01_RawdataStat");
open Out,">$dOut/01_RawdataStat/RawdataStat.xls ";##Sample_ID      Total_Reads     Total_Bases     Error%  Q20%    Q30%    GC%
print Out join("\t","#Sample_ID","Total_Reads","Total_Bases","Q20%","Q30%","GC%"),"\n";
foreach my $sample (sort keys %stat) {
	print Out join("\t",$sample,$stat{$sample}{rawread}/2,$stat{$sample}{rawbase},sprintf("%.2f",100*$stat{$sample}{rawq20}),sprintf("%.2f",100*$stat{$sample}{rawq30}),sprintf("%.2f",100*$stat{$sample}{rawGC})),"\n";
}
close Out;
mkdir "$dOut/02_QC" if (!-d "$dOut/02_QC");
my $dir=basename($qc);
if(!$add){
	open Out,">$dOut/02_QC/QC_stat.xls";
	print Out join("\t","#Sample_ID","Total_Reads","Total_Bases","Q20%","Q30%","GC%"),"\n";
	foreach my $sample (sort keys %stat) {
		print Out join("\t",$sample,$stat{$sample}{cleanread}/2,$stat{$sample}{cleanbase},sprintf("%.2f",100*$stat{$sample}{cleanq20}),sprintf("%.2f",100*$stat{$sample}{cleanq30}),sprintf("%.2f",100*$stat{$sample}{cleanGC})),"\n";
	}
	close Out;
	`ln -s $qc/fig/*.clean*.pdf $dOut/02_QC/`;
#}else{
#	open Out,">$dOut/02_QC/QC_stat.xls";
#	open IN,"$qc/qc"
#	`ln -s $qc/length/*length_distribution* $dOut/02_QC/`;`
}
`ln -s $qc/fig/*.raw*.pdf $dOut/01_RawdataStat/`;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


sub USAGE {
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: caixia.tian

Usage:
  Options:
		-qc	<dir>	fqstp qc dir
		-pollution	<file>	pollution summary
		-o	<dir>	QC report
		-add	if smallRNA
USAGE
	print $usage;
	exit;}
