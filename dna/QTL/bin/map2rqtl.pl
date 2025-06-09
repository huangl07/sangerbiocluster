#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($loc,$map,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"l:s"=>\$loc,
	"m:s"=>\$map,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($loc and $map and $fOut);
open In,$map;
my $groupID;
my %Marker;
my %LG;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^;/) ;
	if (/group/) {
		$groupID=(split(/\s+/,$_))[1];
		$groupID=~s/\D//g;
	}else{
		my ($id,$pos)=split(/\s+/,$_);
		$Marker{$id}=join(",",$groupID,$pos);
	}
}
close In;
open In,$loc;
open Out,">$fOut";
my @head;
my %info;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/);
	if (/MarkerID/) {
		my (undef,$head)=split(/\s+/,$_,2);
		$head=~s/\t/,/g;
		@head=split(/\,/,$head);
		print Out "Genotype,,,",join(",",@head),"\n";
	}else{
		my ($id,$info)=split(/\s+/,$_,2);
		$info=~s/\t/,/g;
		$info=~s/X/H/g;
		$info=~s/U/-/g;
		next if(!exists $Marker{$id});
		print Out join(",",$id,$Marker{$id},$info),"\n"
	}
}
close In;
close Out;
open Out,">$fOut.phe";
print Out join(",","Genotype",@head),"\n";
my @trt;
for(my $i=0;$i<@head;$i++){
	push @trt,1;
}
print Out join(",","id",@trt),"\n";
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"l:s"=>\$loc,
	"m:s"=>\$map,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
