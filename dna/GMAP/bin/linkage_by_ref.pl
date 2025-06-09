#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
$Script=~s/\.pl//g;
my @Times=localtime();
my $year=$Times[5]+1990;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($marker,$fOut,$thread,$chr);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$marker,
				"o:s"=>\$fOut,
				"chr:s"=>\$chr,
				) or &USAGE;
&USAGE unless ( $fOut and $marker );

open In,$marker;
my %chr;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/ ||/MarkerID/);
    my ($id,@info)=split(/\s+/,$_);
    my @a=split(/\-/,$_,2);
    $chr{$a[0]}{$a[0]}=1;
}
close In;
open Out,">$fOut";
foreach my $chr(sort keys %chr){
	print Out ">$chr\n";
	print Out join("\t",keys %{$chr{$chr}}),"\n";
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

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub USAGE {#
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[$month:$day:$year:]
	Contact:Huang Long <huangl\@biomarker.com.cn>
	Options:
				"i:s"=>\$marker,
				"o:s"=>\$fOut,
				"t:s"=>\$thread,
				"c:s"=>\$chr
		-h	Help

USAGE
	print $usage;
	exit;
}
