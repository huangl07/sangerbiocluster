#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($in,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$in,
	"o:s"=>\$out,
			) or &USAGE;
&USAGE unless ($out);
$in||="fastq.gz";
my @file=glob("*.$in");
#if(scalar@file eq "0"){
#	@file=glob("*.fastq");
#}
open Out,">$out";
print Out "MD5\tFile Name\n";
foreach my $file (@file) {
	my $md5sum=`md5sum $file`;
	print Out $md5sum;
}
close Out;

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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
  -i	<format>	input file format(fastq.gz|fastq|vcf)
  -o	<file>	output file name
  -h         Help
USAGE
        print $usage;
        exit;
}
