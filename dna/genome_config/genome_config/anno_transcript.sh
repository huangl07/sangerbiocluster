#!/bin/bash

# @Last-edit Time 2022/12/1
# @Author yiwei.tang
# @mail yiwei.tang@majorbio.com

Usage(){
    echo "Usage:"
    echo "anno transcript.sh -g ref.gff -a anno.summary -o output"
    echo "Description:"
    echo "-a anno.summary, 功能注释"
    echo "-o 输出文件前缀"
    exit -1
}

while getopts ":a:o:" opt_name
do
    case $opt_name in
        a) ANNO="$OPTARG"
            ;;
        o) OUTPUT="$OPTARG"
            ;;
        :) Usage;;
    esac
done

echo -e "chr\ttranscript_start\ttranscript_end\ttranscript_id\tko_id\tko_anno\tgo_term\tgo_anno" > ${OUTPUT}

awk -v FS="\t" -v OFS="\t" 'NR>1{split($1,trans,":"); print trans[3],trans[4],trans[5],trans[2],$6,$7,$8,$9}' $ANNO |sort -u >> ${OUTPUT}
