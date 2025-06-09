#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($cns,$group,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"l:s"=>\$cns,
	"g:s"=>\$group,
	"o:s"=>\$dOut,
			) or &USAGE;
&USAGE unless ($cns and $group and $dOut);
mkdir $dOut if(!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
open In,$cns;
my %cns;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my @a=split(/\s+/,$_);
    $cns{$a[0]}=$a[1];
}
close In;
open In,$group;
my %filehand;
my %head;
open Out,">$dOut/all.segment.txt";
open SH,">$dOut/gistic.sh";
print SH "/mnt/lustre/users/sanger-dev/app/bioinfo/dna/gistics/gistic2 -b ./ -fname all -seg $dOut/all.segment.txt -refgene ~/app/bioinfo/dna/gistics/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat && Rscript ~/app/bioinfo/dna/dna/gistic.R --groupid $dOut/all --out ./\n";
my $head=0;
while(<In>){
    chomp;
    next  if ($_ eq ""||/^$/);
    my ($sample,$gid)=split(/\s+/,$_);
    if(!exists $filehand{$gid}){
        open $filehand{$gid},">$dOut/$gid.segment.txt";
        $head{$gid}=0;
        print SH "/mnt/lustre/users/sanger-dev/app/bioinfo/dna/gistics/gistic2 -b ./ -fname $gid -seg $dOut/$gid.segment.txt -refgene ~/app/bioinfo/dna/gistics/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat && Rscript ~/app/bioinfo/dna/dna/gistic.R --groupid $dOut/$gid --out ./\n";
    }
    open fIn,$cns{$sample};
    my $read;
    while($read=<fIn>){
        chomp $read;
        next if($read eq "" || $read =~ /^$/);
        $read=~s/chr//g;
        my @a=split(/\s+/,$read);
        if($read=~/start/){
            if($head ==0){
                print Out join("\t","sampleID\t$a[0]\t$a[1]\t$a[2]\t$a[6]\t$a[4]\n");
                $head=1;
            }
            if($head{$gid}==0){
                print {$filehand{$gid}} join("\t","sampleID\t$a[0]\t$a[1]\t$a[2]\t$a[6]\t$a[4]\n");
                $head{$gid}=1;
            }
            next;
        }
        print {$filehand{$gid}} join("\t","$sample\t$a[0]\t$a[1]\t$a[2]\t$a[6]\t$a[4]\n");
        print Out join("\t","$sample\t$a[0]\t$a[1]\t$a[2]\t$a[6]\t$a[4]\n");
    }
    close fIn;
}
close Out;
close In;
foreach my $key(sort keys %filehand){
    close $filehand{$key};
}
close SH;
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
	"l:s"=>\$cns,
	"g:s"=>\$group,
	"o:s"=>\$dOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
