#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($dIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"d:s"=>\$dIn,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($dIn and $fOut);
my $pwd=`pwd`;
my $input=`readlink -f $dIn`;
if($pwd eq $input){
	die "you must out the input dir!!!"
}

my $files=`find $dIn/ -type f`;
my @file=split(/\n/,$files);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
open Out,">$fOut";
open LOG,">$fOut.log";
$dIn=`readlink -f $dIn`;
chomp $dIn;
my $dirname=(split(/\//,$dIn))[-2];
my $dir=($year+1900)."-".($mon+1)."-".$mday;
foreach my $file(@file){
    $file=~s/^\.\///g;
	my $filesize=(split(/\s+/,`du $file`))[0];
	print $file,"\t",$filesize,"\n";	
	if($filesize/1024/1024 > 50){
		print "split -b 20G $file $file.";
		`split -b 20G $file $file.`;
		my @file_upload=glob("$file.*");
		foreach my $fs (@file){
			print "aws s3 cp $fs s3://majorbio-public/$dir/$fs\n";
			my $log=`aws s3 cp $fs s3://majorbio-public/$dir/$fs`;
			print LOG $log;
	    	print Out "http://awscdn.majorbio.com/$dir/$fs\n";
		}
	}else{
    	print "aws s3 cp $file s3://majorbio-public/$dir/$file\n";
    	my $log=`aws s3 cp $file s3://majorbio-public/$dir/$file`;
    	print LOG $log;
    	print Out "http://awscdn.majorbio.com/$dir/$file\n";
	}
}
close LOG;
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:$0
Version: $version
Contact: huangl <long.huang\@majorbio.com>
		
Usage:
  Options:
  -d   <dir>    upload dir 
  -o   <file>	download list file
  -h         Help

USAGE
	print $usage;
	exit;
}
