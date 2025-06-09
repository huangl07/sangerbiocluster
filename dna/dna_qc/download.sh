#!/bin/sh
list=$1
dir="download"
n=0
mkdir -p $dir
rm -f $dir/fq.list
cut -f 5,37 $list | sed 1d | sed '/^$/d' | while read id path
do
    echo $id
    echo $path
    n=$((n+1))
    IFS=';' read -r -a array <<< "$path"
    mkdir -p $dir/${n}_${id}
    FileTransfer -type download -s3_path "${array[0]}" -t $dir/${n}_$id
    FileTransfer -type download -s3_path "${array[1]}" -t $dir/${n}_$id
    echo -e "$id\t$dir/${n}_$id/$(basename ${array[0]})\t$dir/${n}_$id/$(basename ${array[1]})" >> $dir/fq.list
done

if [ $2 == "yes" ]; then
    rm -rf final
    mkdir -p final
    echo "Sample Initial Name	Sample Analysis Name	Batch	Library	File Name" > final/fqlist.txt
    cat $dir/fq.list | while read id fq1 fq2
    do
        if [[ ! $(cut -f 1 $dir/fq.list | grep ${id} | wc -l) -gt 1 ]]; then
            cp -l $fq1 final/${id}_R1.fq.gz
            cp -l $fq2 final/${id}_R2.fq.gz
        else
            cat $fq1 >> final/${id}_R1.fq.gz
            cat $fq2 >> final/${id}_R2.fq.gz
        fi
        echo -e "$id\t$id\tbatch1\tWGS\t${id}_R1.fq.gz,${id}_R2.fq.gz" >> final/fqlist.txt.tmp
    done
    sort -u final/fqlist.txt.tmp >> final/fqlist.txt
fi

### 说明
## 用法 bash -e download.sh <样本导出表> yes/no
## 本脚本用于下载fastq文件，输入文件为样本导出表，第二个参数为是否需要合并文件，yes为合并，no为不合并
