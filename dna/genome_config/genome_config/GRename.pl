#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fGenome,$fOut,$Gff,$match);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Copy qw(copy cp);

my $version="1.0.1";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fGenome,
	"o:s"=>\$fOut,
	"g:s"=>\$Gff,
	"f:s"=>\$match
	) or &USAGE;
&USAGE unless ($fGenome and $fOut);
my %change;
if ($match) {
	open In,$match;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$level,$newid)=split(/\s+/,$_);
		$change{$id}=$newid
	}
	close In;

	# 处理fa
	open In,$fGenome;
	open Out,">$fOut.fa";
	$/=">";
	my $fn=0;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ ||/^#/);
		my ($info,@seq)=split(/\n/,$_);
		my $id=(split(/\s+/,$info))[0];
		die "fa && match file doesn't match\n" if (!exists $change{$id});
		print Out ">$change{$id}\n";
		print Out join("\n",@seq),"\n";
	}
	close In;
	close Out;

	if ($Gff) {
		# 处理gff
		open In,$Gff;
		open Out,">$fOut.gtf";
		$/="\n";
		while (<In>) {
			chomp;
			next if ($_ eq ""||/^$/||/^#/);
			my ($id,@info)=split(/\t/,$_);
			die "gff && match file doesn't match\n" if (!exists $change{$id});
			print Out join("\t",$change{$id},@info),"\n";
		}
		close In;
		close Out;
	}
}else{
	copy($fGenome,"$fOut.fa") or die "Copy failed: $!";
	if ($Gff) {
		copy($Gff,"$fOut.gtf") or die "Copy failed: $!";
	}
}



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        yiwei.tang\@majorbio.com;
Script:			$Script
Description:
	reformat genome,rename scaffold name at genome fa file and gff file
eg:
	perl $Script -i Genome.fa -g Genome.gff -k keyname -o dir

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -g	<file>	input genome gff file,
  -o	<str>	output file prefix
  -f	<file>	chromosome change file

  -h         Help

USAGE
        print $usage;
        exit;
}
