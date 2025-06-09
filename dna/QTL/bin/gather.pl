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
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
$fIn=ABSOLUTE_DIR($fIn);

mkdir $fOut if(!-d $fOut);
if(-d "$fOut/result"){
    system("rm -rf $fOut/result");
}
mkdir "$fOut/result" if(!-d "$fOut/result");

my @qtl=glob("$fIn/02.qtl/*/*.detail.result");
my @btl=glob("$fIn/03.btl/*/*.detail.result");
my @enrich=glob("$fIn/06.enrich/*.region_list.txt");

foreach my $qtl(@qtl){
    my $dirname=(split("/",dirname($qtl)))[-1];
    my $basename=basename($qtl);
    my $dir=dirname($qtl);
    mkdir "$fOut/result/$dirname" if(!-d "$fOut/result/$dirname");
    mkdir "$fOut/result/$dirname/qtl/" if(!-d "$fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.png $fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.pdf $fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.result $fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.csv $fOut/result/$dirname/qtl/");
}
foreach my $qtl(@btl){
    my $dirname=(split("/",dirname($qtl)))[-1];
    my $basename=basename($qtl);
    my $dir=dirname($qtl);
    mkdir "$fOut/result/$dirname" if(!-d "$fOut/result/$dirname");
    mkdir "$fOut/result/$dirname/qtl/" if(!-d "$fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.png $fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.pdf $fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.result $fOut/result/$dirname/qtl/");
    system("ln -s $dir/*.csv $fOut/result/$dirname/qtl/");
}
    ##publishDir "${params.outdir}/07.result/${method}/enrich", pattern:"${method}.region.gene.xls"
    ##publishDir "${params.outdir}/07.result/${method}/enrich", pattern:"${method}.region.variant.xls"
    ##publishDir "${params.outdir}/07.result/${method}/enrich/KEGG_result", pattern:"*KEGGenrichment*"
    ##publishDir "${params.outdir}/07.result/${method}/enrich/GO_result", pattern:"*GOenrichment*"
foreach my $enrich(@enrich){
    my $dirname=basename($enrich);
    $dirname=~s/\.region_list\.txt//g;
    my $dir=dirname($dirname);
    mkdir "$fOut/result/$dirname/enrich" if(!-d "$fOut/result/$dirname/enrich");
    my @KEGG=glob("$fIn/06.enrich/$dirname/KEGG_result/*enrichment*");
    if(scalar @KEGG > 0){
        mkdir "$fOut/result/$dirname/enrich/KEGG_result/";
        system("ln -s $fIn/06.enrich/$dirname/KEGG_result/*enrichment*    $fOut/result/$dirname/enrich/KEGG_result/");
    }
    my @GO=glob("$fIn/06.enrich/$dirname/GO_result/*enrichment*");
    if(scalar @GO > 0){
        mkdir "$fOut/result/$dirname/enrich/GO_result/";
        system("ln -s $fIn/06.enrich/$dirname/GO_result/*enrichment*    $fOut/result/$dirname/enrich/GO_result/");
    }
    if(scalar @GO > 0 || scalar @KEGG > 0){
        system("ln -s $fIn/06.enrich/$dirname.region.{gene,variant}.xls    $fOut/result/$dirname/enrich/");
    }
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
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
