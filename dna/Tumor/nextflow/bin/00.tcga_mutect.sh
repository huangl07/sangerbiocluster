#!/bin/bash
#@File    :   00.tcga_mutect.sh
#@Time    :   2023/0407
#@Author  :   xiaoya.ye
#@Version :   1.0
#@Contact :   xiaoya.ye@majorbio.com
Usage(){
    echo "Usage:"
    echo "00.tcga_mutect.sh -v vcf"
    echo "Description:"
    echo "vcf"
    exit -1
}
while getopts ":v:" opt_name
do
    case $opt_name in
        v) vcf="$OPTARG"
            ;;
        :) Usage
            ;;
    esac
done
source  ~/app/bioinfo/dna/new.rc
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%TGT\t]\t[%AD\t]\t%ANN\n' -H ${vcf}| perl -ne 'chomp;if(/#/){print "${_}\n"}else{%stat={};@a=split(/\t/,$_);@b=(split(/\,/,$a[-1]));@c=split(/\,/,$a[3]);foreach $b(@b){$d=(split(/\|/,$b))[0];$b=~s/\|/\t/g;next if(exists $stat{$d});if($c[0] eq $d || $c[1] eq $d ){if(!exists $stat{$d}){$stat{$d}=$b}}};print join("\t",@a[0..$#a-1],join("\t;\t",sort values %stat)),"\n"}' > tcga_mutect.result
exit 0
