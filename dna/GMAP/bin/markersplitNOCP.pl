#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$dOut,$chr,$size);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"fOut:s"=>\$fOut,
    "dOut:s"=>\$dOut,
    "chr:s"=>\$chr,
    "size:s"=>\$size,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $dOut );
if(!-d $dOut){mkdir $dOut}
$dOut=ABSOLUTE_DIR($dOut);
my %chr;
if(defined $chr){
    open In,$chr;
    while(<In>){
        chomp;
        next if ($_ eq ""||/^$/);
        $chr{$_}=1;
    }
    close In;
}
my %chrlen;
if(defined $size){
    open In,$size;
    while(<In>){
        next if ($_ eq ""||/^$/);
        my ($chr,$len)=split(/\s+/,$_);
        $chrlen{$chr}=$len;
    }
    close In;
}
open In,$fIn;
my @Indi;
my $head;
my %filehand;
open Out,">$fOut";
my %sca;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    if(/^#MarkerID/ ||/MarkerID/){
        (undef,undef,@Indi)=split(/\s+/,$_);
        $head=join("\t","chr","pos",@Indi);
    }else{
        s/aa/a/g;
        s/bb/b/g;
        s/ab/h/g;
        s/--/-/g;
        my ($id,$type,@indi)=split(/\s+/,$_);
        my ($chr,$pos)=split(/\-/,$id);
        if (!exists $sca{$chr} && !exists $chr{$chr}){
            $sca{$chr}=1;
        }
        if(!exists $filehand{$chr}){
            open $filehand{$chr},">$dOut/$chr.marker\n";
            print Out $chr,"\t",int($chrlen{$chr}/(scalar @indi)/2/20),"\t$dOut/$chr.marker","\n";
            print {$filehand{$chr}} $head,"\n";
        }
            print {$filehand{$chr}} join("\t",$chr,$pos,join("\t",@indi)),"\n"; 
        
    }
}
close In;
close Out;
open Out,">$dOut/compare.sh";
if(scalar keys %sca != 0){
    my @chr=sort keys %filehand;
    for(my $i=0;$i<@chr;$i++){
        for(my $j=$i+1;$j<@chr;$j++){
            my $sca=$chr[$i];
            my $chr=$chr[$j];
            next if($sca eq $chr);
            next if (scalar keys %chr !=0  && exists $chr{$sca} && exists $chr{$chr});
            next if(scalar keys %chr !=0 && !exists $chr{$sca} && !exists $chr{$chr});
            print Out "-i1 $sca.bin.marker -i2 $chr.bin.marker -o $sca-$chr.mlod\n";
        }
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
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -input	<file>	input file name
  -fout	<file>	input list file
  -dout <dir>   input dir file
  -chr  <file>  chr file for compare
  -h         Help

USAGE
        print $usage;
        exit;
}
