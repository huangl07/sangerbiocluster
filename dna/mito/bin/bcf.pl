#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bcf,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"b:s"=>\$bcf,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($bcf);
open In,"bcftools query -f \"%CHROM\t%POS\t%REF\t%ALT\t[%TGT\t]%DP4\n\" $bcf|";
my @fasta="";
open Out,">$fOut.variant.table";
print Out join("\t","CRHOM","POS","REF","ALT","$fOut.GT","$fOut.AD"),"\n";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my @a=split(/\t/,$_);
    my @b=split(/\,/,$a[-1]);
    my @allele=split(",",join(",",$a[2],$a[3]));
    my @gt=split(/\//,$a[4]);
    my $ad1=$b[0]+$b[1];
    my $ad2=$b[2]+$b[3];
    if($gt[0] ne $gt[1] || $gt[0] ne $a[2]){
        print Out join("\t",@a[0..$#a-1],"$ad1,$ad2"),"\n";
    }
    if($ad1 >= $ad2){
        push @fasta,$gt[0];
    }else{
        push @fasta,$gt[1];
    }
}
close In;
close Out;
open Out,">$fOut.consensus.fa";
print Out ">$fOut\_chrM\n",join("",@fasta),"\n";
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
