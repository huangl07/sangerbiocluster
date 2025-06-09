#!/usr/bin/perl -w
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
	"in:s"=>\$in,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($in and $out);
######################################################
open IN,$in;
open OUT,">$out";
my %hash;
while(<IN>){
	chomp;
	next if ($_ eq ""||/^$/);
    my ($sample,$gro) = split/\t/,$_;
    if (!defined $gro){
		print "can't get groupname\n";
	}
	push @{$hash{$gro}},$sample;
}
close IN;
foreach (sort keys %hash){
    print OUT "$_:\t",join("\t",@{$hash{$_}}),"\n";
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR
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
Contact:	tong.wang\@majorbio.com
Version:	$version
Script:		$Script
Description:	make new group.list for vcf2treemix.py
Usage:
  Options:
  -in	<file>	input raw group file 
  -out	<file>	output new group file 
  -h		Help

USAGE
        print $usage;
        exit;
}
