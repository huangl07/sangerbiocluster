#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$map);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
        "input:s"=>\$fIn,
        "m:s"=>\$map,
	"output:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$map;
my %info;
my $lgID;
while (<In>) {
        chomp;
	next if ($_ eq ""||/^$/ || /^;/ );
        if(/group/){
                $lgID=(split(/\s+/,$_))[-1];
                next;
        }
	my ($id,$pos)=split(/\s+/,$_);
	$info{$id}{lgID}=$lgID;
	$info{$id}{pos}=$pos;
}
close In;
open In,$fIn;
open Out,">$fOut";
open Draw,">$fOut.draw";

while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    s/kh/hk/g;
    s/pn/np/g;
    s/ml/lm/g;
    my ($marker,$type,$phase,@geno)=split(/\s+/,$_);
    next if(!exists $info{$marker});
    $phase=~s/{|}//g;
    $type=~s/<|>//g;
    $phase=~s/1/2/g;
    $phase=~s/0/1/g;

    my @out;
    my @out1;
    my @out2;
    if(!defined $info{$marker}{lgID}){print $marker;die;}
    if(!defined $info{$marker}{pos}){print $marker;die;}
    push @out1,join(",",$marker,$info{$marker}{lgID},$info{$marker}{pos});
    push @out2,join(",",$marker,$info{$marker}{lgID},$info{$marker}{pos});
    for(my $i=0;$i<@geno;$i++){
 
        my $haplosource=determineHaploSource($type,$phase,$geno[$i]);
        my @hap=split("",$haplosource);
        if(scalar @hap  != 2){print "$geno[$i]\t$haplosource\n","\n";die;}
        push @out,$haplosource;
        push @out1,$hap[0];
        push @out2,$hap[1];
    }
    print Out join("\t",$marker,"<$type>","{$phase}",@out),"\n";
    print Draw join(",",@out1),"\n";
    print Draw join(",",@out2),"\n";
}
close In;
close Out;
close Draw;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub determineHaploSource {#
        my ($type,$linkPhase,$progenyGeno)=@_;

        return "--" if ($progenyGeno eq "--") ;

        my $haploSource='';;
        my (%parentAllel,@gameteCombination)=();

        my %haploIndex=(
                "0"=>"11",
                "1"=>"12",
                "2"=>"21",
                "3"=>"22",
        );

        my ($p1,$p2,$m1,$m2)= $type =~/(.)(.)x(.)(.)/ ;
        @{$parentAllel{"P"}}=($p1,$p2);
        @{$parentAllel{"M"}}=($m1,$m2);

        my ($PLinkPhase, $MLinkPhase) = split //,$linkPhase ;
        if ($PLinkPhase eq '1'){
                @{$parentAllel{"P"}} = reverse @{$parentAllel{"P"}};
        }
        if ($MLinkPhase eq '1'){
                @{$parentAllel{"M"}} = reverse @{$parentAllel{"M"}} ;
        }

        foreach my $pGamete (@{$parentAllel{"P"}}) {
                foreach my $mGamete (@{$parentAllel{"M"}}) {
                        push @gameteCombination,$pGamete.$mGamete;
                }
        }

        for (my $j=0;$j<@gameteCombination ;$j++) {
                
                my @haplo=split //,$gameteCombination[$j];
                my $haplo=join("",sort {$a cmp $b} @haplo);

                if ($haplo eq $progenyGeno) {

                        my ($allel1,$allel2) = $progenyGeno =~/(.)(.)/;

                        if ($p1 eq $p2) {

                                my ($lp) = $haploIndex{$j} =~/.(.)/;
                                $haploSource = "0".$lp;
                                last;
                                
                        }elsif ($m1 eq $m2) {

                                my ($lp) = $haploIndex{$j} =~/(.)./;
                                $haploSource = $lp."0";
                                last;

                        }elsif (join("",sort {$a cmp $b} ($p1,$p2)) eq join("",sort {$a cmp $b} ($m1,$m2)) and $allel1 ne $allel2) {
                                
                                $haploSource = "--";
                                last;

                        }else{
                                
                                $haploSource = $haploIndex{$j};
                        }
                }
                
        }
        return $haploSource;

}


sub haplo2gene {
	my ($cross_type,$phase,@haplo)=@_;
	my $gene;
	my @gene;
	$cross_type=~s/>|x|<//g;
	my @cross_type = split //,$cross_type;
	my @phase = split //,$phase;
	for (my $i=0;$i<2;$i++) {
		if ($phase[$i] eq '1') {
			my $t=$cross_type[$i*2];
			$cross_type[$i*2] = $cross_type[$i*2+1];
			$cross_type[$i*2+1] = $t;
		}
		if ($haplo[$i] eq '1') {
			$gene[$i]=$cross_type[$i*2];
		}elsif($haplo[$i] eq '2'){
			$gene[$i]=$cross_type[$i*2+1];
		}else{
			if ($cross_type[$i*2] eq $cross_type[$i*2+1]) {
				$gene[$i] = $cross_type[$i*2];
			}else{
				$gene[$i]='-';
				$gene='--';
			}
		}
	}
	my @unit= @gene;
	if (!defined $gene) {
		$gene=$unit[0].$unit[1];
	}
	return $gene;
}
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
    "input:s"=>\$csv,
	"output:s"=>\$fOut,
  -h         Help

USAGE
        print $usage;
        exit;
}
