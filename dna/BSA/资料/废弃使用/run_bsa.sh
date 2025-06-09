#!bin/bash
###
 # @Author: error: git config user.name && git config user.email & please set dead value or install git
 # @Date: 2022-12-06 13:17:21
 # @LastEditors: error: git config user.name && git config user.email & please set dead value or install git
 # @LastEditTime: 2023-01-16 17:39:55
 # @FilePath: /ZQH_5612/mnt/ilustre/users/jing.wang/work/FX2022112400238/test/run_BSA.sh
 # @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
### 
Usage(){
    echo "Usage:"
    echo "bash run_BSA.sh -t type"
    echo "Description:"
    echo "type: resume,要不要重跑"
    exit -1
}

while getopts ":t:" opt_name
do
    case $opt_name in
        w) winsize="$OPTARG"
            ;;
        s) stepsize="$OPTARG"
            ;;
        p) pvalue="$OPTARG"
            ;;
        t) type="$OPTARG"
            ;;
        :) Usage;;

    esac
done

if [ ! $winsize ];then
    winsize=1000000
fi

if [ ! $stepsize ];then
    stepsize=10000
fi

if [ ! $pvalue ];then
    pvalue=0.001
fi


if [ ! -d ./result ];then
    mkdir result
fi

if [ ! -d ./data ];then
    mkdir data
fi


workdir=`pwd`

##根据不同情况进行调整参数
if [ ! $type ];then
    nextflow run -bg /mnt/ilustre/users/yuan.xu/scripts/pipeline/BSA_production/BSA/BSA.nf --vcf $workdir/workflow_results/05.annovar/combine_variants/pop.final.vcf --outdir $workdir/result --chr $workdir/data/chr.list --gff $workdir/workflow_results/02.reference/ref.gtf --anno_summary $workdir/workflow_results/02.reference/anno.summary --group $workdir/data/group.list --winsize $winsize --stepsize $stepsize --pvalue $pvalue
else
    nextflow run -bg /mnt/ilustre/users/yuan.xu/scripts/pipeline/BSA_production/BSA/BSA.nf --vcf $workdir/workflow_results/05.annovar/combine_variants/pop.final.vcf --outdir $workdir/result --chr $workdir/data/chr.list --gff $workdir/workflow_results/02.reference/ref.gtf --anno_summary $workdir/workflow_results/02.reference/anno.summary --group $workdir/data/group.list --winsize $winsize --stepsize $stepsize --pvalue $pvalue -resume
fi
