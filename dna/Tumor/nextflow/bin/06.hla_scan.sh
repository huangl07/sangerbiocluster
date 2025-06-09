#!/bin/bash
#@File    :   06.hypermutation.sh
#@Time    :   2023/04/21
#@Author  :   xiaoya.ye
#@Version :   1.0
#@Contact :   xiaoya.ye@majorbio.com
set +e
Usage(){
    echo "Usage:"
    echo "06.hypermutation.sh -b sample_bam -s sample"
    echo "Description:"
    echo "hla_scan"
    exit -1
}
while getopts ":b:s:l:" opt_name
do
    case $opt_name in
        b) sample_bam="$OPTARG"
            ;;
        s) sample="$OPTARG"
            ;;
        l) hlalist="$OPTARG"
            ;;
        :) Usage
            ;;
    esac
done
source  ~/app/bioinfo/dna/new.rc
touch ${sample}.report 
xargs -a ${hlalist} -I {} -P 4 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hla_dataset/hla_scan_r_v2.1.4 -v 38 -d /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hla_dataset/db/HLA-ALL.IMGT -g {} -b ${sample_bam} -t 8 -f 75 >> ${sample}.report || exit 0