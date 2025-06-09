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
my ($mlod,$fOut,$thread,$chr);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$mlod,
				"o:s"=>\$fOut,
				"t:s"=>\$thread,
				"c:s"=>\$chr
				) or &USAGE;
&USAGE unless ($mlod and $fOut and $thread );
open In,$chr;
my %chr;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	my @a=split;
	push @{$chr{$a[0]}},$a[0];
}
close In;
my %LG;
open In,$mlod;
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/ || /MLOD/);
	my ($m1,$m2,$lod)=split(/\,/,$_);
	my $sca=$m1;
	my $chr=$m2;
	if(exists $chr{$sca}){
		$sca=$m2;
		$chr=$m1;
	}
	if(!exists $LG{$sca} ||$LG{$sca}{max}<$lod){
		$LG{$sca}{max}=$lod;
		$LG{$sca}{chr}=$chr;
		$LG{$sca}{num}=1;
	}elsif($LG{$sca}{max} == $lod){
		$LG{$sca}{max}=$lod;
		$LG{$sca}{chr}=$chr.",".$LG{$sca}{chr};
		$LG{$sca}{num}++;
	}
}
close In;
open Out,">$fOut.stat";
foreach my $sca(sort keys %LG){
	if($LG{$sca}{num} == 1 && ${LG{$sca}{max} > $thread}){
		push @{$chr{$LG{$sca}{chr}}},$sca;
		print Out join("\t",$sca,$LG{$sca}{max},$LG{$sca}{num},$LG{$sca}{chr},"passed"),"\n";
	}else{
		print Out join("\t",$sca,$LG{$sca}{max},$LG{$sca}{num},$LG{$sca}{chr},"failed"),"\n";
	}
}
close Out;
open Out,">$fOut";
foreach my $chr(sort keys %chr){
	print Out ">$chr\n";
	print Out join("\t",@{$chr{$chr}}),"\n";
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
				"i:s"=>\$mlod,
				"o:s"=>\$fOut,
				"t:s"=>\$thread,
				"c:s"=>\$chr
		-h	Help

USAGE
	print $usage;
	exit;
}
