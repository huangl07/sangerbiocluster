#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dIn,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"d:s"=>\$dIn,
	"o:s"=>\$dOut,
			) or &USAGE;
&USAGE unless ($dIn and $dOut);
my @bams=glob("$dIn/*.depth");
my %stat;
foreach my $bam(@bams){
    my $sample=basename($bam);
    $sample=~s/\.depth//g;
    open In, "$bam";
    my $n=0;
    while(<In>){
        chomp;
        next if ($_ eq ""||/^$/||/^#/);
        #rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
        my @a=split;
        $stat{$sample}{tags}++ if($a[3] !=0);
        $stat{$sample}{meandepth}+=$a[6];
    }
    close In;
}
my @file=glob("$dIn/populations.log.distribs");
open In, $file[0];
my $flag="NULL";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/ ||/^#/ );
    if(/END/){
        $flag="NULL";
        next;
    }
    if(/BEGIN loci_per_sample/){
        $flag="loci_per_sample";
        <In>;
        next;
    }
    if(/BEGIN variant_sites_per_sample/){
        $flag="variant_sites_per_sample";
        <In>;
        next;
    }
    if($flag eq "NULL"){
        next;
    }
    if($flag eq "loci_per_sample"){
        #sample  n_loci  present_loci    missing_loci    frequency_missing
        next if (/n_loci/);
        my @a=split;
        $stat{$a[0]}{ntags}=$a[2];
        $stat{$a[0]}{totaltags}=$a[1];
    }
    if($flag eq "variant_sites_per_sample"){
        #sample  n_sites present_sites   missing_sites   frequency_missing
        next if (/n_loci/);
        my @a=split;
        $stat{$a[0]}{nvariant}=$a[2];
        $stat{$a[0]}{totalvariant}=$a[1];
    }
}
close In;
open Out,">$dOut/sample.stat";
print Out "#sampleID\ttotal_tags\ttotal_snps\tpresent_tags\tpresent_snps\tmean_depth\n";
foreach my $sample(sort keys %stat){
    next if ($sample eq "sample");
    $stat{$sample}{totaltags}||=0;
    $stat{$sample}{totalvariants}||=0;
    $stat{$sample}{nvariant}||=0;
    $stat{$sample}{ntags}||=0;
    $stat{$sample}{meandepth}||=0;
    if($stat{$sample}{meandepth} == 0){
        print Out join("\t",$sample,$stat{$sample}{totaltags}, $stat{$sample}{totalvariant},$stat{$sample}{ntags},$stat{$sample}{nvariant},sprintf("%.2f",0)),"\n";
    }else{
        print Out join("\t",$sample,$stat{$sample}{totaltags}, $stat{$sample}{totalvariant},$stat{$sample}{ntags},$stat{$sample}{nvariant},sprintf("%.2f",$stat{$sample}{meandepth}/$stat{$sample}{tags})),"\n";
    }
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
	"i:s"=>\$dIn,
	"o:s"=>\$dOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
