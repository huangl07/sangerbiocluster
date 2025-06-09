#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$dfile,$compare);
use Data::Dumper;
use JSON;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "d:s"=>\$dfile,
	"o:s"=>\$fOut,
	"c:s"=>\$compare,
			) or &USAGE;
&USAGE unless ($dfile and $fOut);
my @file=glob("$dfile/*.vcf.gz");
my %group;
my %chr;
foreach my $vcf(@file){
    my $filename=basename($vcf);
    $filename=~s/.vcf.gz//g;
    my ($chr,$gid)=split(/\-/,$filename);
    $group{$gid}{$chr}=ABSOLUTE_DIR($vcf);
    $chr{$chr}=1;
}
my @chr=sort keys %chr;
open In,$compare;
open Out,">$fOut";
while(<In>){
	chomp;
	next if ($_ eq ""||/^$/);
	my ($g1,$g2)=split(/\t/,$_);
	for(my $c=0;$c<@chr;$c++){
		print Out "$group{$g1}{$chr[$c]}\t$group{$g2}{$chr[$c]}\t$g1\_$g2\t$chr[$c]\n";
	}
}
close In;
close Out;
#open Out,">$fOut";
#my @gid=sort keys %group;
#my @chr=sort keys %chr;
#for (my $i=0;$i<@gid;$i++){
#    for(my $j=$i+1;$j<@gid;$j++){
#        for(my $c=0;$c<@chr;$c++){
#            print Out "$group{$gid[$i]}{$chr[$c]}\t$group{$gid[$j]}{$chr[$c]}\t$gid[$i]\_$gid[$j]\t$chr[$c]\n";
#        }
#    }
#}
#close Out;

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
