#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"result:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out,">$fOut.region";
open Out1,">$fOut.xls";
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/ );
    if(/lod/){my @a=split;print Out1 join("\t",$a[-1],@a[0..$#a-1]),"\n";next;}
    #lod.chr lod.pos lod.ci.low lod.ci.high lod.lod pos1 pos2 var qname
#2 54.4805381158822 18.0208151316231 56.2641907988019 5.60748343487518 chr1-38880198_39792549 chr1-13106480_15258292 2.86642610623377 blup
#7 19.0703472361186 15.9902123959584 21.4775213820959 7.25476639144608 chr6-103036125_116400643 chr6-68749917_98843620 7.32925113528473 blup
#10 43.2769942742431 40.3906538137461 44.9214240157752 8.65020363643208 chr9-127150802_129667388 chr9-131743500_145816815 5.7564700653841 blup
    my ($chr,$peek,$pos1,$pos2,$lod,$mark1,$mark2,$var,$qname)=split(/\s+/,$_);
    my @pos;
    if($mark1 =~ /([^\-]*)\-(\d+)\_(\d+)/){
        push @pos,$2;
        push @pos,$3;
    }else{
        my (undef,$pos)=split(/\-/,$mark1);
        push @pos,$pos;
    }
    if($mark2 =~ /([^\-]*)\-(\d+)\_(\d+)/){
        push @pos,$2;
        push @pos,$3;
    }else{
        my (undef,$pos)=split(/\-/,$mark2);
        push @pos,$pos;
    }

    @pos=sort{$a<=>$b} @pos;
    print Out join("\t",$chr,$pos[0],$pos[-1],$qname),"\n";
    print Out1 join("\t",$qname,$chr,$peek,$pos1,$pos2,$lod,$pos[0],$pos[-1],$var),"\n";
}
close In;
close Out;
close Out1;
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
	"result:s"=>\$fIn,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
