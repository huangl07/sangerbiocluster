#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$dmap,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($dmap and $fOut);
my @maps=glob("$dmap/*.result.map");
	my %Marker;

foreach my $map(@maps){
	open In,$map;
	my $gid;
	$gid=(split(/\./,basename($map)))[0];
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ ||/^;/) ;
		if (/group/) {
			next;
		}else{
			my ($id,$female,$male,$sexAver)=split(/\s+/,$_);
			if ($female ne "NA") {
				$Marker{$gid}{female}{$id}=$female;
			}
			if ($male ne "NA") {
				$Marker{$gid}{male}{$id}=$male;
			}
			if ($sexAver ne "NA") {
				$Marker{$gid}{sexAver}{$id}=$sexAver;
			}
		}
	}
	close In;
}

mkdir  $fOut if(!-d $fOut);
open SexAver,">$fOut/total.sexAver.map";
open Male,">$fOut/total.male.map";
open Female,">$fOut/total.female.map";
foreach my $group (keys %Marker){
	print SexAver "group $group\n";
	print Male "group $group\n";
	print Female "group $group\n";
	foreach my $id(sort {$Marker{$group}{female}{$a} <=>$Marker{$group}{female}{$b}} keys %{$Marker{$group}{female}}){
		print Female $id,"\t",$Marker{$group}{female}{$id},"\n";
	}
	foreach my $id(sort {$Marker{$group}{male}{$a} <=>$Marker{$group}{male}{$b}} keys %{$Marker{$group}{male}}){
		print Male $id,"\t",$Marker{$group}{male}{$id},"\n";
	}
	foreach my $id(sort {$Marker{$group}{sexAver}{$a} <=>$Marker{$group}{sexAver}{$b}} keys %{$Marker{$group}{sexAver}}){
		print SexAver $id,"\t",$Marker{$group}{sexAver}{$id},"\n";
	}
}
close Male;
close Female;
close SexAver;
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
  -i	<file>	input cross link file 
  -o	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}
