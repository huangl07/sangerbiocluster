#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Scalar::Util qw(looks_like_number);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
    "popt:s"=>\$popt
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
if($popt eq "CP" || $popt eq "F1"){
    open In,$fIn;
    my %trt;
    my @trt;
    my $nind=0;
    while (<In>){
        chomp;
        next if ($_ eq ""||/^$/);
        if(/sampleID/){
            (undef,@trt)=split(/\s+/,$_);
            for(my $i=0;$i<@trt;$i++){
               push @{$trt{$trt[$i]}},join(" ","sampleID",$trt[$i])
            }
            next;
        }
        $nind++;
        my @line=split(/\s+/,$_);
        for(my $i=1;$i<@line;$i++){
            if(!looks_like_number($line[$i])){$line[$i] = "*"}
            if($line[$i] eq "NA"){$line[$i]="*"}
            push @{$trt{$trt[$i-1]}},join(" ",$line[0],$line[$i])
        }
    }
    close In;
    open List,">$fOut.list";
    foreach my $trt(@trt){
        open Out,">$fOut.$trt.txt";
        print Out "nind=$nind\nntrt=2\nmiss=*\n";
		#print Out "name $trt\n";
        print Out join("\n",@{$trt{$trt}}),"\n";
        close Out;
        print List $trt,"\t",ABSOLUTE_DIR("$fOut.$trt.txt"),"\n";
    }
    close List;
}else{
    open In,$fIn;
    my %trt;
    my @trt;
    my $nind=0;
    my @indi;
    while (<In>){
        chomp;
        next if ($_ eq ""||/^$/);
        if(/sampleID/){
            (undef,@trt)=split(/\s+/,$_);
            next;
        }
        my @line=split(/\s+/,$_);
        push @indi,$line[0];
        $nind++;
        for(my $i=1;$i<@line;$i++){
            if(!looks_like_number($line[$i])){$line[$i] = "*"}
			if($line[$i] eq "NA"){$line[$i]="*"}
            push @{$trt{$trt[$i-1]}},$line[$i];
        }
    }
    close In;
    open List,">$fOut.list";
    foreach my $trt(@trt){
        open Out,">$fOut.$trt.txt";
        print Out join(",","Genotype",@indi),"\n";
        print Out join(",",$trt,@{$trt{$trt}}),"\n";
        close Out;
        print List $trt,"\t",ABSOLUTE_DIR("$fOut.$trt.txt"),"\n";
    }
    close List;
}
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
