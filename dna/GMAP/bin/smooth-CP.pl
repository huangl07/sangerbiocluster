#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.1";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($map,$loc,$fKey,$dOut,$win,$ami);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$map,
				"l:s"=>\$loc,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				"win:s"=>\$win,
				"ami:s"=>\$ami,
				) or &USAGE;
&USAGE unless ($map and $loc and $fKey and $dOut);
#------------------------------------------------------------------
# Global parameter  
#------------------------------------------------------------------
$dOut||="./";
$win||=30;
$ami||=1;
mkdir $dOut unless (-d $dOut) ;
open In,$loc;
my %pri;
my $nind;
my $head;
my %primary;
my $nloc=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^;/);
	if (/\=/ || /^#/) {
		$head.="$_\n";
	}else{
		my ($id,$type,$phase,@geno)=split(/\s+/,$_);
		$nloc=scalar @geno;
		$type=~s/<|>//g;
		$phase=~s/\{|\}//g;
		$nind=scalar @geno;
		$phase=~s/1/2/g;
		$phase=~s/0/1/g;
		$pri{$id}{phase}=$phase;
		$pri{$id}{type}=$type;

		for (my $i=0;$i<@geno;$i++) {
			my $haplosource=determineHaploSource($type,$phase,$geno[$i]);
			if ($haplosource eq "") {
				$haplosource = "--";
			}
			push @{$pri{$id}{geno}},$haplosource;
			push @{$primary{$id}},$geno[$i];
		}
	}
}
close In;
open In,$map;
my %order;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /group/);
	my ($id,$pos)=split(/\s+/,$_);
	next if (!exists $pri{$id});
	$order{$id}=$pos;
}
close In;

my @order=sort {$order{$a}<=> $order{$b}} keys %order;
my $step=1;
my %correct;
for (my $i=0;$i<@order;$i++) {
	my $start=$i-$win;
	my $end=$i+$win;
	$start=0 if ($start < 0);
	$end=scalar @order -1 if ($end > scalar @order -1);
	for (my $j=0;$j<@{$pri{$order[$i]}{geno}};$j++) {
		my %stat;
		for (my $k=$start;$k<=$end;$k++) {
			if (!exists $pri{$order[$k]}{geno}) {
				print $k;die;
			}
			if (length($pri{$order[$k]}{geno}[$j]) < 2) {
				print $pri{$order[$k]}{geno}[$j],"\t",$order[$k],"\t","$j";die;
			}
			my ($b1,$b2)=split(//,$pri{$order[$k]}{geno}[$j]);
			$stat{1}{$b1}++ if($b1 !~ /-/ && $b1 != 0);
			$stat{2}{$b2}++ if($b2 !~ /-/ && $b2 != 0);
		}
		my @b1=sort{$stat{1}{$b}<=>$stat{1}{$a}} keys %{$stat{1}};
		my @b2=sort{$stat{2}{$b}<=>$stat{2}{$a}} keys %{$stat{2}};
		my ($b1,$b2)=split(//,$pri{$order[$i]}{geno}[$j]);
		if (scalar @b1 > 1 && $stat{1}{$b1[0]}==$stat{1}{$b1[1]}) {
			$b1="0";
		}elsif (scalar @b1 !=0) {
			$b1=$b1[0];
		}elsif (scalar @b2 == 0) {
			$b1="0";
		}

		if (scalar @b2 > 1 && $stat{2}{$b2[0]}==$stat{2}{$b2[1]}) {
			$b2="0";
		}elsif (scalar @b2 !=0) {
			$b2=$b2[0];
		}elsif (scalar @b2 == 0) {
			$b2="0";
		}
		$correct{$order[$i]}{geno}[$j]=join("",$b1,$b2);
	}
}
open Out,">$dOut/$fKey.correct.loc";
open Cor,">$dOut/$fKey.correct.log";
open Phase,">$dOut/$fKey.correct.phase";
#print Out $head,"\n";
my %amibigous;
foreach my $id (sort {$order{$a}<=> $order{$b}} keys %order) {
	my @out;
	my @outl;
	my @outp;
	my $cor=0;
	for (my $i=0;$i<@{$correct{$id}{geno}};$i++) {
		my @hap=split(//,$correct{$id}{geno}[$i]);
		my $genes=haplo2gene($pri{$id}{type},$pri{$id}{phase},$hap[0],$hap[1]);
		if($pri{$id}{type} eq "hkxhk" && $genes eq "--"){
			$genes="hk";
		}
		push @out,$genes;
		if ($pri{$id}{geno}[$i] eq $correct{$id}{geno}[$i] || ($correct{$id}{geno}[$i] eq "--" && $pri{$id}{geno}[$i] ne "--" )) {
			push @outp,$pri{$id}{geno}[$i];
		}else{
			push @outp,$pri{$id}{geno}[$i]."->".$correct{$id}{geno}[$i];
		}
		if($primary{$id}[$i] ne $genes){
			$cor++;
		}
		push @outl,$correct{$id}{geno}[$i];
	}
	$pri{$id}{phase}=~s/1/0/g;
	$pri{$id}{phase}=~s/2/1/g;
	$amibigous{$id}{stat}=$cor/$nloc;
	if($cor/$nloc > $ami){
		$amibigous{$id}{type}="T";
		print Out join("\t",$id,"<".$pri{$id}{type}.">","{".$pri{$id}{phase}."}",@{$primary{$id}}),"\n";
		print Cor join("\t","$id(amibiguous)","<".$pri{$id}{type}.">","{".$pri{$id}{phase}."}",@{$pri{$id}{geno}}),"\n";
		print Phase join("\t",$id,"<".$pri{$id}{type}.">","{".$pri{$id}{phase}."}",@outl),"\n";
	}else{
		print Out join("\t",$id,"<".$pri{$id}{type}.">","{".$pri{$id}{phase}."}",@out),"\n";
		print Cor join("\t",$id,"<".$pri{$id}{type}.">","{".$pri{$id}{phase}."}",@outp),"\n";
		print Phase join("\t",$id,"<".$pri{$id}{type}.">","{".$pri{$id}{phase}."}",@outl),"\n";
	}

	
}
close Out;
close Cor;
close Phase;
open Out,">amibiguous.log";
print Out "Total: ",scalar keys %pri,"\n";
foreach my $id(sort keys %amibigous){
	$amibigous{$id}{type}||="F";
	print Out join("\t",$id,$amibigous{$id}{stat},$amibigous{$id}{type}),"\n";
}
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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
	$cross_type=~s/x//g;
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
	my @unit=sort {$a cmp $b} @gene;
	if (!defined $gene) {
		$gene=$unit[0].$unit[1];
	}
	return $gene;
}


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version	
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription: 
	-window_type����Խ��-diff_ratio��ֵԽ�ͣ�����Խ�ֱ���
	
Usage:
  Options:
  -m		<file>	Map file, forced
  -l		<file>	Loc file with linkage phase info,forced

  -k		<str>	Key of output file,forced
  -d		<str>	dir of output file,default ./
  -win		<num>	win size 30
  -h		Help

USAGE
	print $usage;
	exit;
}
