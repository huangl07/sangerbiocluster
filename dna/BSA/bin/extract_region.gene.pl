#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw ($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($region,$gff,$table,$fOut,$extype);
GetOptions(
    "help|?"=>\&USAGE,
    "region:s"=>\$region,
    "table:s"=>\$table,
    "out:s"=>\$fOut,
    "type:s"=>\$extype,
) or &USAGE;
&USAGE unless ($region and  $table and $fOut);
$extype||="gene";
my %region;
open IN,$region;
while(<IN>){
    chomp;
    next if ($_ eq "" ||/^$/ );
	#    next if (/chr/);
    my ($chr,$pos1,$pos2,undef)=(split/\s+/,$_,4);
	next if ($chr eq "chr");
    $region{$chr}{join("\t",$pos1,$pos2)}=1;
}
close IN;
my %stat;
open Out,">$fOut.total.variant";
open IN,$table;
my %variant;
while(<IN>){
    chomp;
    if (/CHROM/){
        next;
    }else{
        my @lines=split(/\s+/,$_);
        my $ANN=$lines[-1];
        my $chr=$lines[0];
        my $pos=$lines[1];
        my $type=$lines[4];
        my ($impact,$genes)=(split(/\|/,$ANN))[2,4];
        foreach my $region (sort keys %{$region{$chr}}){
            my ($pos1,$pos2)=split(/\t/,$region);
            if ($pos>=$pos1 && $pos<=$pos2){
                
                if ($type eq "SNP"){
                    $stat{$chr}{$region}{snp}++;
                    #$stat{$chr}{$region}{effgene}{$genes} =1  if($impact eq "HIGH" || $impact eq "MODERATE");
                    #$stat{$chr}{$region}{effsnp}++ if($impact eq "HIGH" || $impact eq "MODERATE");
                    push @{$variant{$chr}{$region}},"$chr,$region\t".$_;
                }
                if ($type eq "INDEL"){
                    $stat{$chr}{$region}{indel}++;
                    #$stat{$chr}{$region}{effgene}{$genes} =1 if($impact eq "HIGH" || $impact eq "MODERATE");
                    #$stat{$chr}{$region}{effindel}++ if($impact eq "HIGH" || $impact eq "MODERATE");
                    push @{$variant{$chr}{$region}},"$chr,$region\t".$_;

                }
            }
        }

    }
}
close IN;
foreach my $chr(keys %variant){
    foreach my $region(sort keys %{$variant{$chr}}){
        print Out join("\n",@{$variant{$chr}{$region}}),"\n";
    }
}
close Out;
# open IN,$gff;
#     while(<IN>){
#         chomp;
#         next if (/^$/ || /^\#/);
#         my ($chr,$type,$start,$end)=(split(/\t/,$_))[0,2,3,4];
#         next if ($type ne "mRNA");
#         foreach my $region (sort keys %{$region{$chr}}){
#             my ($pos1,$pos2)=split(/\s+/,$region);
#             if ($end>$pos1 && $start<$pos2){
#                 $stat{$chr}{$region}{gene}++;
#             }
#         }
#     }  

# close IN;
# open OUT,">$fOut.total.variant.stat";
# print OUT "ChromosomeID\tStart\tEnd\t#SNP\t#Indel\t#Effsnp\t#EffInDel\n";
# foreach my $chr (sort keys %region){
#     foreach my $region(sort keys %{$region{$chr}}){
#         $stat{$chr}{$region}{snp}||=0;
#         $stat{$chr}{$region}{effsnp}||=0;
#         $stat{$chr}{$region}{indel}||=0;
#         $stat{$chr}{$region}{effindel}||=0;
#         print OUT join ("\t",$chr,$region,$stat{$chr}{$region}{snp},$stat{$chr}{$region}{indel},$stat{$chr}{$region}{effsnp},$stat{$chr}{$region}{effindel}),"\n";
#     }
# }
# close OUT;
# open OUT,">$fOut.gene.stat";
# print OUT "ChromosomeID\tStart\tEnd\t#Gene\t#GeneEff\n";
# foreach my $chr (sort keys %region){
#     foreach my $region(sort keys %{$region{$chr}}){
#         $stat{$chr}{$region}{gene}||=0;
#         my $num=0;
#         if(!exists $stat{$chr}{$region}{effgene}){
#             $num=0;
#         }else{
#             $num = scalar keys %{$stat{$chr}{$region}{effgene}};
#         }
#         print OUT join("\t",$chr,$region,$stat{$chr}{$region}{gene},$num ),"\n";
#     } 
# }
# close OUT;
########################################################################
print STDOUT "\nDone. Total elasped time :",time()-$BEGIN_TIME,"s\n";
########################################################################
sub USAGE{
    my $usage=<<"USAGE";
Contact:        dong.wang\@majorbio.com;
Script:         $Script
Description:
    eg:
    perl $Script -gff -region -table -out
Usage:
    Options:
    -gff          <file>        gff3 file
    -region       <file>        region file
    -table        <file>        variant file
    -out          <file>        output prefix name
    -h             Help          
USAGE
    print $usage;
    exit;
}
