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
my ($fqlist,$dOut,$dShell,$queue,$majorbio);
GetOptions(
			"help|?" =>\&USAGE,
			"fqlist:s"=>\$fqlist,
			"out:s"=>\$dOut,
			"dsh:s"=>\$dShell,
			"queue:s"=>\$queue,
			"majorbio:s"=>\$majorbio,
				) or &USAGE;
&USAGE unless ($fqlist and $dOut and $dShell);
#$queue||="dna";
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
mkdir "$dOut/qc" if(!-d "$dOut/qc");
$fqlist=ABSOLUTE_DIR($fqlist);
$dOut=ABSOLUTE_DIR($dOut);
#my $qsub="qsub-slurm.pl";
#my $fastp="fastp";
open In,$fqlist;
open SH,">$dShell/step01.fastp-trim.sh";
mkdir "$dOut/fig" if (!-d "$dOut/fig");
open Out,">$dOut/fastq.list";
open Stat,">$dOut/stat.list";
open Json,">$dOut/json.list";
my %stat;
my $num;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$fq1,$fq2,undef)=split(/\s+/,$_);
	$sample=~s/-/\_/g;
	push @{$stat{$sample}{fq1}},$fq1;
	push @{$stat{$sample}{fq2}},$fq2;
}
close In;

foreach my$sample(sort keys%stat){
	if(scalar@{$stat{$sample}{fq1}} ne "1"){
		print SH "cat ",join("\t",@{$stat{$sample}{fq1}})," >$dOut/$sample.raw.1.fastq.gz && ";
		print SH "cat ",join("\t",@{$stat{$sample}{fq2}})," >$dOut/$sample.raw.2.fastq.gz && ";
	}else{
		print SH "mv  ${$stat{$sample}{fq1}}[0] $dOut/$sample.raw.1.fastq.gz && ";
		print SH "mv  ${$stat{$sample}{fq2}}[0] $dOut/$sample.raw.2.fastq.gz && ";
	}
	print SH "fastp  -i $dOut/$sample.raw.1.fastq.gz -I $dOut/$sample.raw.2.fastq.gz -o $dOut/$sample.clean.1.fastq.gz -O $dOut/$sample.clean.2.fastq.gz ";
	print SH "--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  --adapter_sequence_r2  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT " if($majorbio == 1);
	print SH " -q 20 --length_required 30 -5 20 -3 20 -W 4 -M 20 -n 10 -j $dOut/qc/$sample.json -h $dOut/qc/$sample.html -w 8 && ";

	print SH " perl $Bin/bin/fastp.pl -i $dOut/qc/$sample.json -o $dOut/qc/$sample && Rscript $Bin/bin/ngsqc.r --base $dOut/qc/$sample.raw.atgcn --qual $dOut/qc/$sample.raw.qual --key $sample.raw --od $dOut/fig && Rscript $Bin/bin/ngsqc.r --base $dOut/qc/$sample.clean.atgcn --qual $dOut/qc/$sample.clean.qual --key $sample.clean --od $dOut/fig \n";

	print Out "$sample\t$dOut/$sample.clean.1.fastq.gz\t$dOut/$sample.clean.2.fastq.gz\n";
	print Stat "$sample\t$dOut/qc/$sample.stat\n";
	print Json "$sample\t$dOut/qc/$sample.json\n";
}
close In;
close SH;
close Out;
close Stat;
close Json;
my $job="qsub-slurm.pl $dShell/step01.fastp-trim.sh --Resource mem=12G --CPU 8 --Maxjob 40 --Queue $queue";
`$job`;

open SH2,">$dShell/step01.2.md5sum.sh";
print SH2 "cd $dOut && perl $Bin/bin/md5sum.pl -i raw.*.fastq.gz -o raw.md5.txt \n";
print SH2 "cd $dOut && perl $Bin/bin/md5sum.pl -i clean.*.fastq.gz -o clean.md5.txt \n";
close SH2;
`qsub-slurm.pl --Resource mem=8G --CPU 1 $dShell/step01.2.md5sum.sh --Queue $queue`;
print "md5sum done!\n";

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

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: caixia.tian

Usage:
  Options:
	-fqlist	<file>	input file 
	-out	<dir>	output dir
	-dsh	<dir>	output work sh
	-queue	<str>	DNA or RNA ,defult RNA
	-majorbio <num> 1 for majorbio library ;2 for others
USAGE
	print $usage;
	exit;
}
