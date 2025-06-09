#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fconfig,$dOut,$fqlist,$step,$stop,$runID,$hostname,$queue,$project_info,$majorbio);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"outdir:s"=>\$dOut,
	"fqlist:s"=>\$fqlist,
        "project_info:s"=>\$project_info,
		"majorbio:s"=>\$majorbio,
			) or &USAGE;
&USAGE unless ($fqlist);
#$queue||="RNA";
$queue="none";
$hostname = `hostname`;
$queue = "DNA" if ($hostname =~ /^compute-5/);
$queue = "SANGERDEV" if ($hostname =~ /^login/);
if ($queue eq "none") {
    die "Error: Can't detect the queue name, please check the cluster.\n";
}
$majorbio||=1;
$dOut||=`pwd`;
$dOut=(split(/\s+/,$dOut))[0];
mkdir $dOut if (!-d $dOut);
$fqlist=ABSOLUTE_DIR($fqlist);
$project_info=ABSOLUTE_DIR($project_info);
$dOut=ABSOLUTE_DIR($dOut);
mkdir "$dOut/work_sh" if (!-d "$dOut/work_sh");
my $tmp=time();
open Log,">$dOut/work_sh/QC.$tmp.log";
	print Log "########################################\n";
	print Log "fastq trim  !"; my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.fastp-qc.pl -fqlist $fqlist -out $dOut/01.CleanData/ -dsh $dOut/work_sh --queue $queue ";
	if($majorbio == 1){
		$job.= "--majorbio 1"
	}
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log " fastq trim  Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	print Log "merge result!"; $time=time();
	print Log "########################################\n";
	my $qc=ABSOLUTE_DIR("$dOut/01.CleanData/");
	$job="perl $Bin/bin/step02.merge_report.pl -qc $qc -o $dOut/02.QCreport/ ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "merge result Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	print Log "generate qc report" ,time()-$time,"s\n";
	print Log "########################################\n";
	$job="python $Bin/bin/qc_report/Report.py  --project $project_info --result  $dOut ";
	print Log "$job\n";
	`$job`;
	print Log "########################################\n";
	print Log "data release" ,time()-$time,"s\n";
	my $job="sh $Bin/bin/release.sh  $dOut  ";
	print Log "$job\n";
	`$job`;
	print Log "########################################\n";
	print Log "merge result Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
close Log;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
  -fqlist	<dir>	input fastq list
  -outdir	<dir>	output file dir
  -queue	<str>	DNA or SANGERDEV
  -majorbio	<num>	1 majorbio library and 2 others
  -h         Help

USAGE
        print $usage;
        exit;
}
