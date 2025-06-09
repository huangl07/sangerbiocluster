#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
$Script=~s/\.pl//g;
my @Times=localtime();
my $year=$Times[5]+1990;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$fPosi);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out1,">$fOut.parent1.loc";
open Out2,">$fOut.parent2.loc";
open Out3,">$fOut.phase";
my $head=0;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($id,$type,$phase,@geno)=split(/\s+/,$_);
    if($head == 0){
        my @out;
        for(my $i=0;$i<@geno;$i++){
            push @out,"S$i";
        }
        print Out1 join("\t","chr","pos",@out),"\n";
        print Out2 join("\t","chr","pos",@out),"\n";
        $head=1;
    }
    my ($chr,$pos)=split(/_/,$id);
    my @out1;
    my @out2;
    my @out3;
    my %basesource=("0"=>"-","1"=>"a","2"=>"b","-"=>"-",);
    for (my $i=0;$i<@geno;$i++){
        my $haplosource=determineHaploSource($type,$phase,$geno[$i]);
        my @base=split(//,$haplosource);
        if($type =~ /lmxll/ || $type =~ /hkxhk/){
            push @out1,$basesource{$base[0]};
        }
        if($type =~ /nnxnp/||$type =~ /hkxhk/){
            push @out2,$basesource{$base[1]};
        }
        push @out3,$haplosource;
    }

    if($type =~ /lmxll/ || $type =~ /hkxhk/){
        print Out1 $chr,"\t",$pos,"\t",join("\t",@out1),"\n";
    }
     if($type =~ /nnxnp/ || $type =~ /hkxhk/){
        print Out2 $chr,"\t$pos\t",join("\t",@out2),"\n";
    }
    print Out3 $chr,"\t$pos\t$type\t$phase\t",join("\t",@out3),"\n";

}
close In;
close Out1;
close Out2;
close Out3;
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
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath{
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chop $return;
                $return .="/".$file;
                chdir($pwd);
        }
        return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[$month:$day:$year:]
	Contact:Huang Long <huangl\@biomarker.com.cn>
	Options:
		-i	<dir>	input file
		-o	<file>	output file
		-h	Help

USAGE
	print $usage;
	exit;
}
