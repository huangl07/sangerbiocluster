#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gff , $anno , $fOut ,$select);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$select,
	"a:s"=>\$anno,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($select and $anno and $fOut );
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:

Usage:
  Options:
  -i	<file>  input select file name
  -a	<file>  input anno file name
  -o	<key>   output keys of file name
  -h         Help

USAGE
        print $usage;
        exit;
}
my $number||=1;
open In,$select;
my %region;
my $qtl;
while (<In>) {
	chomp;
	s/\"//g;
	next if ($_ eq ""||/^$/ ||/#/ ||/mark1/);
	next if ($. == 1);
	my ($chr,$mark1,$mark2)=split(/\s+/,$_);
#	my $start1=$mark1;
#	my $start2=$mark2;
#	$region{$chr}{$number}=join("\t",sort{$a<=>$b}($start1,$start2));
	$region{$chr}{$number}=join("\t",sort{$a<=>$b}($mark1,$mark2));
	$number++ ;
}
close In;

open In,$anno;
my %edetail;
my %info;
my %einfo;
my %stat;
my %astat;
my $head;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/);
	if (/^#/) {
		$head=$_;
		next;
	}#Gene_name	Gene_id	Transcript_id	Bio_Type	Chr	Pos1	Pos2	High	Moderate	Low	Modifier	NR-ID	NR-ANNO	Uni-ID	Uni-ANNO	KEGG-ID	KEGG-ANNO	GO-ID	GO-ANNO	EggNOG-ID	EggNOG-ANNO
	my ($Gene_name,$Gene_id,$Transcript_id,$Bio_Type,$chr,$Pos1,$Pos2,$High,$Moderate,$Low,$Modifier,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$express)=split(/\t/,$_);
	#	next if($chr !~ /chr/);
	#	$chr=(split(/\D+/,$chr))[-1];
	my $regioned=0;
	next if(!exists $region{$chr});
	foreach my $number (sort{$a<=>$b} keys %{$region{$chr}}){
		my ($pos3,$pos4)=split(/\t/,$region{$chr}{$number});
		if (($Pos1 > $pos3 && $Pos1 < $pos4) || ($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
				push @{$info{$chr}{$number}},$_;
				$stat{$chr}{$number}{qtl}=$region{$chr}{$number};
				$stat{$chr}{$number}{total}++;
				$stat{$chr}{$number}{totalnr}++ if($nrid ne "--");
				$stat{$chr}{$number}{totaluni}++ if($uniid ne "--");
				$stat{$chr}{$number}{totalkegg}++ if($koid ne "--" && $koid ne "");
				$stat{$chr}{$number}{totalgo}++ if($goid ne "--" && $goid ne "");
				$stat{$chr}{$number}{totaleggnog}++ if($egid ne "S");
				$stat{$chr}{$number}{eff}++ if($High+$Moderate > 0);
				push @{$einfo{$chr}{$number}},$_ if($High+$Moderate > 0);
				$stat{$chr}{$number}{effnr}++ if($nrid ne "--" && $High+$Moderate > 0);
				$stat{$chr}{$number}{effuni}++ if($uniid ne "--" && $High+$Moderate > 0 );
				$stat{$chr}{$number}{effkegg}++ if( ($koid ne "--" && $koid ne "") && $High+$Moderate > 0);
				$stat{$chr}{$number}{effgo}++ if(($goid ne "--" && $goid ne "") && $High+$Moderate > 0);
				$stat{$chr}{$number}{effeggnog}++ if($egid ne "S" && $High+$Moderate > 0);
				$regioned=1;
		}
	}
	$astat{total}++;
	$astat{eff}++ if($regioned == 1);
}
close In;
open Out,">$fOut.total";
print Out "#\@chr\tpos1\tpos2\ttotal\teff\ttotalnr\ttotaluni\ttotalkegg\ttotalgo\ttotaleggnog\teffnr\teffuni\teffkegg\teffgo\teffeggnog\n";
print Out "$head\n";
foreach my $chr (sort{$a cmp $b} keys %region) {
	foreach my $number (sort{$a<=>$b} keys %{$stat{$chr}}) {
		$stat{$chr}{$number}{totalnr}||=0;
		$stat{$chr}{$number}{totaluni}||=0;
		$stat{$chr}{$number}{totalkegg}||=0;
		$stat{$chr}{$number}{totalgo}||=0;
		$stat{$chr}{$number}{totaleggnog}||=0;
		$stat{$chr}{$number}{effnr}||=0;
		$stat{$chr}{$number}{effuni}||=0;
		$stat{$chr}{$number}{effkegg}||=0;
		$stat{$chr}{$number}{effgo}||=0;
		$stat{$chr}{$number}{effeggnog}||=0;
		$stat{$chr}{$number}{total}||=0;
		$stat{$chr}{$number}{eff}||=0;
		print Out join("\t","\@$chr",$stat{$chr}{$number}{qtl},$stat{$chr}{$number}{total},$stat{$chr}{$number}{eff},$stat{$chr}{$number}{totalnr},$stat{$chr}{$number}{totaluni},$stat{$chr}{$number}{totalkegg},$stat{$chr}{$number}{totalgo},$stat{$chr}{$number}{totaleggnog},$stat{$chr}{$number}{effnr},$stat{$chr}{$number}{effuni},$stat{$chr}{$number}{effkegg},$stat{$chr}{$number}{effgo},$stat{$chr}{$number}{effeggnog}),"\n";
		next if (!defined $info{$chr}{$number});
		print Out join("\n",@{$info{$chr}{$number}}),"\n";
	}
}
close Out;
open Out,">$fOut.eff";
print Out "#\@chr\tpos1\tpos2\n";
print Out "$head\n";
foreach my $chr (sort{$a cmp $b} keys %region) {
	foreach my $number (sort keys %{$region{$chr}}) {
		print Out "\@$chr\t$region{$chr}{$number}\n";
		next if (!defined $einfo{$chr}{$number});
		print Out join("\n",@{$einfo{$chr}{$number}}),"\n";
#		foreach my $line (@{$info{$chr}{$number}}) {
#			my @split=split(/\t/,$line);
#			print Out $line,"\n" if($split[7]+$split[8] > 0);
#		}
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
