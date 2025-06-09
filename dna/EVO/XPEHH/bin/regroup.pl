#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$group,$map);
use Data::Dumper;
use JSON;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
    "g:s"=>\$group,
    "m:s"=>\$map,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$map;
open Out,">$fOut.map";
while(<In>){
    chomp;
    next if ($_ eq ""||/^S/);
    my @a=split;
    $a[0]=~s/\D//g;
    print Out join(" ",@a),"\n";
}
close Out;
close In;
open In,$group;
my %group;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my @a=split(/\s+/,$_);
    $group{$a[0]}=$a[1];
}
close In;
open In,$fIn;
open Out,">$fOut.ped";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my @a=split;
    print Out join(" ",$group{$a[0]},@a[1..$#a]),"\n";
}
close Out;
close In;
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
  -i	<file>	input dict name
  -o	<file>	out qc base file name
  -h         Help

USAGE
        print $usage;
        exit;
}
