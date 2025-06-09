#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$popt,$PID,$MID);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"PID:s"=>\$PID,
    "MID:s"=>\$MID,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out,">$fOut";
my $Ppos=-1;
my $Mpos=-1;
my %fstat;
while(<In>){
    chomp;
    next if ($_ eq ""||/^$/);
    if(/^#/){
        my ($id,@indi)=split(/\s+/,$_);
        my @out;
        push @out,"#MarkerID\ttype";
        for (my $i=0;$i< @indi;$i++){
			if($indi[$i]=~/\]([^:]*)\:/){
                $id=$1;
                print $id,"\n";
                if($id eq $PID){
                    $Ppos=$i;
                }elsif($id eq $MID){
                    $Mpos=$i;
                }else{
                    push @out,$id;
                }
            }else{
                $id=$indi[$i];
                if($id eq $PID){
                    $Ppos=$i;
                }elsif($id eq $MID){
                    $Mpos=$i;
                }else{
                    push @out,$id;
                }
            }
        }
        die "$PID\t$Ppos\t$MID\t$Mpos\n" if($Ppos < 0 || $Mpos < 0);
        print Out join("\t",@out),"\n";
    }else{
        my ($id,@indi)=split(/\s+/,$_);
        if($indi[$Ppos] eq "./." || $indi[$Mpos] eq "./." || $indi[$Ppos] eq "--" || $indi[$Mpos] eq "--"){
            $fstat{PMmissing}++;
            next;
        }
        my ($p1,$p2)=split("/",$indi[$Ppos]);
        my ($m1,$m2)=split("/",$indi[$Mpos]);
        my %stat;
        $stat{$p1}++;
        $stat{$p2}++;
        $stat{$m1}++;
        $stat{$m2}++;
        my %geno;
        my $type;
        if(scalar keys %stat == 1 && $p1 eq $p2 && $m1 eq $m2 && $p1 eq $m1){#aaxaa
            $type="aaxaa";
            #$geno{$p1}="a";
            $fstat{aaxaa}++;
            next;
        }elsif(scalar keys %stat == 2 && $p1 eq $p2 && $m1 eq $m2 ){#aaxbb
            $type="aaxbb";
            $geno{$p1}="a";
            $geno{$m1}="b";
        }elsif(scalar keys %stat == 2 && $p1 eq $p2 && $m1 ne $m2 ){#nnxnp
            $type="nnxnp";
            my ($g1,$g2)=sort {$stat{$b} <=> $stat{$a}} keys  %stat;
            $geno{$g1}="n";
            $geno{$g2}="p"
        }elsif(scalar keys %stat ==2 && $p1 ne $p2 && $m1 eq $m2){#lmxll
            $type="lmxll";
            my ($g1,$g2)=sort {$stat{$b} <=> $stat{$a}} keys  %stat;
            $geno{$g1}="l";
            $geno{$g2}="m"
        }elsif($p1 ne $p2 && $m1 ne $m2 && $p1 eq $m1 && $p2 eq $m2){#hkxhk
            $type="hkxhk";
            $geno{$p1}="h";
            $geno{$p2}="k";
        }elsif(scalar keys %stat == 3 && $p1 ne $p2 && $m1 ne $m2){#efxeg
            $type="efxeg";
            my ($g1,$g2,$g3)=sort {$stat{$b} <=> $stat{$a}} keys  %stat;
            $geno{$g1}="e";
            if($g1 eq $p1){
                $geno{$p2} = "f";
            }else{
                $geno{$p1} = "f";
            }
            if($g1 eq $m1){
                $geno{$m2} = "g";
            }else{
                $geno{$m1} = "g";
            }
        }elsif(scalar keys %stat ==4 ){#abxcd
            $type="abxcd";
            $geno{$p1}="a";
            $geno{$p2}="b";
            $geno{$m1}="c";
            $geno{$m2}="d"
        }elsif(scalar keys %stat == 3 && $p1 ne $p2 && $m1 eq $m2 ){#abxcc
            $type="abxcc";
            $geno{$p1}="a";
            $geno{$p2}="b";
            $geno{$m1}="c"
        }elsif(scalar keys %stat == 3 && $p1 eq $p2 && $m1 ne $m2 ){#ccxab
            $type="ccxab";
            $geno{$p1}="c";
            $geno{$m1}="a";
            $geno{$m2}="b"
        }else{
            print $p1,"\t",$p2,"\t",$m1,"\t",$m2;
            print $_;die;
            next;
        }
        my @out;
        push @out,$id;
        push @out,$type;
        $fstat{$type}++;
        for(my $i=0;$i<@indi;$i++){
            next if ($i == $Ppos || $i == $Mpos);
            my ($g1,$g2)=split("/",$indi[$i]);
            if($indi[$Ppos] eq "./." || $indi[$Ppos] eq "--" || !exists $geno{$g1} || !exists $geno{$g2}){
                push @out,"--";
            }else{
                my $gtype=join("",sort($geno{$g1},$geno{$g2}));
                if($gtype eq "pp"){$gtype = "np"};
                if($gtype eq "mm"){$gtype = "lm"};
                push @out,$gtype;
            }
        }
        print Out join("\t",@out),"\n";
    }
}
close Out;
close In;
open Out,">$fOut.stat";
foreach my $type(keys %fstat){
    print Out $type,"\t",$fstat{$type},"\n";
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"PID:s"=>\$PID,
    "MID:s"=>\$MID,
  -h         Help

USAGE
        print $usage;
        exit;
}
