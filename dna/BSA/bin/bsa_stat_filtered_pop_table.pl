#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$fChr,$fanno);
GetOptions(
                "help|?" =>\&USAGE,
                "infile:s"=>\$fIn,
                "outfile:s"=>\$fOut,
                "chrfile:s"=>\$fChr,
                "anno:s"=>\$fanno
                ) or &USAGE;
&USAGE unless ($fIn and $fOut and $fChr);
print STDOUT "正在读取染色体列表...","\n";
my %stat; # 抓取染色体列表
my @chrlist;
open In,$fChr;
while(<In>){
    chomp;
    next if(/^$/ || $_ eq "");
    my ($ochr,undef)=split;
    push @chrlist, $ochr;
    $stat{$ochr}{"snp"}=0;
    $stat{$ochr}{"esnp"}=0;
    $stat{$ochr}{"indel"}=0;
    $stat{$ochr}{"eindel"}=0;
};
close In;
print STDOUT "正在读取anno.summary以提取基因列表...","\n";
my (%genes,$gene_id); #抓取基因列表
open In, $fanno;
readline In;
while(<In>){
    chomp;
    next if(/^$/ || $_ eq "");
    ($gene_id, undef) = split(/:/);
    $genes{$gene_id} = 0;
}
close In;
open In, $fIn;
my $header;
chomp($header = readline In);
my (undef,undef,undef,undef,undef,@Header) = split(/\t/,$header,-1);
print STDOUT "正在读取pop.table完成统计...","\n";
print STDOUT join("\t",@Header),"\n";
my $annIdx=0; # 抓取ANN所在列
my $err=1;
for(my $i=0;$i<@Header;$i++){
    if ($Header[$i] eq "ANN"){
        $annIdx=$i;
        $err=0;
        last;
    };
};
if($err > 0 ){
    warn "No ANN col found!";
    exit;
}
my ($chr, $pos, $ref, $alt, $vtype, @fields,$impact,$geneid,@afields);
while(<In>){
    chomp;
    next if(/^$/ || $_ eq "");
    ($chr, $pos, $ref, $alt, $vtype, @fields) = split;
    next if(!exists $stat{$chr});
    if($vtype eq "SNP"){
        $stat{$chr}{"snp"}++;
    }elsif($vtype eq "INDEL"){
        $stat{$chr}{"indel"}++;
    }else{
        next;
    }
    if ($fields[$annIdx] =~ /\|HIGH\||\|MODERATE\|/){
        if($vtype eq "SNP"){
            $stat{$chr}{"esnp"}++;
        }else{
            $stat{$chr}{"eindel"}++;
        }
    }
    my @ann = split(/;/,$fields[$annIdx],-1);
    for (my $anni=0; $anni<@ann;$anni++){
        @afields = split(/\|/,$ann[$anni]);
        $impact=$afields[2];
        $geneid=$afields[4];
        if (exists $genes{$geneid}){
            $stat{$chr}{"genelist"}{$geneid}++;
            if($impact =~ /HIGH|MODERATE/){
                $stat{$chr}{"egenelist"}{$geneid}++;
            }
        }
    }
}
close In;
open Out, ">$fOut";
print Out join("\t", "CHROM", "SNP_Number", "Effective_SNP", "INDEL_Number", "Effective_INDEL", "Gene_Number", "Effective_Gene"),"\n";
foreach my $c(@chrlist){
    print Out join("\t",
        $c,
        $stat{$c}{"snp"},
        $stat{$c}{"esnp"},
        $stat{$c}{"indel"},
        $stat{$c}{"eindel"},
        scalar keys %{$stat{$c}{"genelist"}},
        scalar keys %{$stat{$c}{"egenelist"}}
    ),"\n";
};
close Out;
print STDOUT "\nDone. Total elasped time :",time()-$BEGIN_TIME,"s\n";
########################################################################
sub USAGE{
    my $usage=<<"USAGE";
Contact:        yiwei.tang\@majorbio.com;
Script:         $Script
Description:
    eg:
    perl $Script --infile <file> --outfile <file> --chrfile <file> --anno <file>
Usage:
    Options:
    -infile          <file>        pop.table file
    -outfile         <file>        output file
    -chrfile         <file>        target chromosome list
    -anno            <file>        anno.summary
    -h             Help
USAGE
    print $usage;
    exit;
}
