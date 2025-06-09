#!/usr/bin/perl

use strict;
use warnings;

my $snpeff;
my @line;
my $DP;
my $AD;
my $minDP;
my $minAF;
my $maf;

my $Hugo_Symbol;
my $Entrez_Gene_Id;
my $Center;
my $NCBI_Build;
my $Chromosome;
my $Start_Position;
my $End_Position;
my $Strand;
my $Variant_Classification;
my $Variant_Type;
my $Reference_Allele;
my $Tumor_Seq_Allele1;
my $Tumor_Seq_Allele2;
my $Tumor_Sample_Barcode;
my $Protein_Change;
my $Matched_Norm_Sample_Barcode;
my $dbSNP_RS;
my $dbSNP_Val_Status;
$snpeff = shift;
$minDP = shift;
$minAF = shift;

($Entrez_Gene_Id,$Center,$NCBI_Build,$Strand)=("NA","NA","hg38","NA");
$Tumor_Sample_Barcode = (split /\./,$snpeff)[0];
$Matched_Norm_Sample_Barcode = (split /\./,$snpeff)[1];
$maf = $Tumor_Sample_Barcode . "_" . $Matched_Norm_Sample_Barcode . ".maf";


open INFILE,"$snpeff" or die "$snpeff";
open OUTFILE,">$maf" or die "$maf";

#print OUTFILE "$Hugo_Symbol\t$Entrez_Gene_Id\t$Center\t$NCBI_Build\t$Chromosome\t$Start_Position\t$End_Position\t$Strand\t$Variant_Classification\t$Variant_Type\t$Reference_Allele\t$Tumor_Seq_Allele1\t$Tumor_Seq_Allele2\t$dbSNP_RS\t$dbSNP_Val_Status\t$Tumor_Sample_Barcode\t$Matched_Norm_Sample_Barcode\n";
my %change_words;
while (<DATA>) {
    chomp;
    my @data = split ("\t", $_);
    $change_words{$data[0]} = $data[1];
}

while(<INFILE>){
    next if /^#/;
    chomp;
    @line = split /\t/;
    next if $line[4] eq '.';
    $DP = (split /:/,$line[-1])[-2];
    next if $DP < $minDP;
    $AD = (split /,/,$line[-1])[-1];
    $dbSNP_RS="";
    $dbSNP_Val_Status="";
    
    ($Chromosome,$Start_Position,$Reference_Allele,$Tumor_Seq_Allele2) = @line[0,1,3,4];
    #$Chromosome =~ s/chr//;
    if ((length($Reference_Allele) == 1) && (length($Tumor_Seq_Allele2) == 1)){
        $End_Position = $Start_Position;
        $Variant_Type ="SNP";
    }elsif ((length($Reference_Allele) >1) && (length($Tumor_Seq_Allele2) == 1)){
        $Start_Position += 1;
        $Reference_Allele = substr($Reference_Allele,1);
        $Tumor_Seq_Allele2 = '-';
        $End_Position = $Start_Position + length($Reference_Allele) - 1;
        $Variant_Type ="DEL";
    }elsif ((length($Reference_Allele) ==1 ) && (length($Tumor_Seq_Allele2) > 1)){
        $End_Position = $Start_Position + 1;
        $Reference_Allele = '-';
        $Tumor_Seq_Allele2 = substr($Tumor_Seq_Allele2,1);
        $Variant_Type ="INS";
    };
    $Tumor_Seq_Allele1 = $Reference_Allele;

    if ($line[7] =~ /(ANN=\S+)/){
        ($Variant_Classification, $Hugo_Symbol)=(split (/\|/, $1))[1,3];
        my @Variant_Classification = split (/&/, $Variant_Classification);
        if (exists $change_words{$Variant_Classification[0]}) {
            $Variant_Classification = $change_words{$Variant_Classification[0]};
        } else {
            print "error: $Variant_Classification not exist\n";
            next;
        }
        next if not $Variant_Classification;
    }else{
        next;
    };
        
    print OUTFILE "$Hugo_Symbol\t$Entrez_Gene_Id\t$Center\t$NCBI_Build\t$Chromosome\t$Start_Position\t$End_Position\t$Strand\t$Variant_Classification\t$Variant_Type\t$Reference_Allele\t$Tumor_Seq_Allele1\t$Tumor_Seq_Allele2\t$dbSNP_RS\t$dbSNP_Val_Status\t$Tumor_Sample_Barcode\t$Matched_Norm_Sample_Barcode\n";
    
};

close INFILE;
close OUTFILE;
__DATA__
missense_variant	Missense_Mutation
synonymous_variant	Nonsense_Mutation
start_retained_variant	Nonsense_Mutation
stop_retained_variant	Nonsense_Mutation
minus_1_frameshift_variant	Frame_Shift_Del
minus_2_frameshift_variant	Frame_Shift_Del
frameshift_variant	Frame_Shift_Ins
inframe_deletion	In_Frame_Del
inframe_insertion	In_Frame_Ins
disruptive_inframe_deletion	In_Frame_Del
disruptive_inframe_insertion	In_Frame_Ins
silent_mutation	Silent
splice_acceptor_variant	Splice_Site
splice_donor_variant	Splice_Site
splice_region_variant	Splice_Site
3_prime_UTR_variant	3'UTR
3_prime_UTR_truncation+exon_loss	3'Flank
5_prime_UTR_variant	5'UTR
5_prime_UTR_truncation+exon_loss_variant	5'Flank
intergenic_region	IGR
conserved_intergenic_variant	IGR
intron_variant	Intron
miRNA	RNA