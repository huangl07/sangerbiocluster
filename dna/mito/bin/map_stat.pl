#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($indir,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"b:s"=>\$indir,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($indir);
my @qc=glob("$indir/*.qc.stat");
open Out,">$fOut";
print Out join("\t","sampleID","rawread","rawbase","rawQ20","rawQ30","rawGC","cleanread","cleanbase","cleanQ20","cleanQ30","cleanGC","capture","chrM-reads","chrM-coverage","chrM-meandepth","allmappedratio","propermappedratio","cov>meaddepth*20%","cov>100","cov>500"),"\n";
foreach my $qc(sort @qc){
	my $id=(split(/\./,basename($qc)))[0];
	my @out;
	push @out,$id;
	open In,$qc;
	my $cleanreads;
	while(<In>){
		chomp;
		next if ($_ eq ""||/^$/ ||/rawreads/);
		my @line=split;
		$cleanreads=$line[7];
		push @out,join("\t",@line[1..5],@line[7..11]);
	}
	close In;

	open In,"$indir/$id.chrM.coverage";
	my $meandepth;
	while(<In>){
		chomp;
		next if($_ eq ""||/^$/ ||/endpos/);
		my @line=split(/\t/,$_);
		push @out,$line[4]/$cleanreads;
		$meandepth=$line[8];
		push @out,join("\t",$line[4],$line[6],$line[7]);
	}
	close In;
	open In,"$indir/$id.depth";
	my %coverage;
	my $sum=0;
	while(<In>){
		chomp;
		next if($_ eq ""||/^$/);
		my ($chr,$pos,$base)=split(/\s+/,$_);
		$sum+=$base;
		$coverage{$meandepth/5}+=$base if($base > $meandepth/5);
		$coverage{100}+=$base if($base >= 100);
		$coverage{500}+=$base if($base >= 500);
	}
	close In;
	open In,"$indir/$id.flagstat";
	my $total;
	my $mapped;
	my $proper;
	while(<In>){
		chomp;
		next if($_ eq ""||/^$/);
		s/%//g;
		if(/0 mapped/){
			my @a=split;
			$mapped=$a[4];
			$mapped =~ s/\(//g;
		}
		if(/properly paired/){
			my @a=split;
			$proper=$a[5];	
			$proper =~ s/\(//g;
		}
	}
	close In;
	$coverage{$meandepth/5}||=0;
	$coverage{100}||=0;
	$coverage{500}||=0;
	$coverage{$meandepth/5}||=0;
	push @out,$mapped;
	push @out,$proper;
	push @out,$coverage{$meandepth/5}/$sum;
	push @out,$coverage{100}/$sum;
	push @out,$coverage{500}/$sum;
	print Out join("\t",@out),"\n"
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"b:s"=>\$indir,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
