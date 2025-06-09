## 插入位点检测SOP

### 需要投递WGS

```
1.	在生信云平台投递WGS任务：
注意，如果没有特别需求，默认都不分析！
2.	投递后，运行结束获得WGS结果路径，必须是完整路径
3.	在172上，工作路径git clone git@git.majorbio.com:long.huang/dna.git 或者
git clone https://git.majorbio.com/long.huang/dna.git
4.	然后git checkout yxy（切换分支）
5.	运行run_tdna.sh
source  ~/app/bioinfo/dna/new.rc
bash $workdir/dna/TDNA/run_tdna.sh $(wgs_v4结果路径) ${ref.fa完整路径} ${insert.fa完整路径}
例如：
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna.sh /mnt/ilustre/isanger_workspaceWgsV4/20230411/WgsV4_719u_6ltlqmhn11iicb8nmabos2 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/ref.fa /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/insert.txt
或者不git，直接使用172有的版本
source  ~/app/bioinfo/dna/new.rc
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna.sh $(wgs_v4结果路径) ${ref.fa完整路径} ${insert.fa完整路径}
例如：
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna.sh /mnt/ilustre/isanger_workspaceWgsV4/20230411/WgsV4_719u_6ltlqmhn11iicb8nmabos2 /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/ref.fa /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/insert.txt
在运行nextflow过程中若出现任何问题，请及时咨询叶晓雅。
如需重新运行意外中断的nextflow，在如上代码后面添加“-resume”即可。
6.	运行tdna报告脚本
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna_report.sh yes/no（是否有WGS运行结果路径）${tdna结果路径} ${WGS的结果路径} 
例如：
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna_report.sh yes /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/test/tang/raw_nf/output /mnt/lustre/sanger_workspaceWgsV4/20230811/WgsV4_03el_gagu2ld1lhk2i24p66r7am
```

### 不需要投递WGS
```
运行nextflow
注：fq.list格式为三列，“\t”分割,分别为raw_data_R1_path, raw_data_R2_path,name
step1： 若git，则使用git clone git@git.majorbio.com:long.huang/dna.git 或者
git clone https://git.majorbio.com/long.huang/dna.git
step2：
使用dna/TDNA路径下的tdna.nf和run_tdna_report.sh即可
***注意！！！
运行nextflow给的outdir需要和git下来的dna文件夹在同一路径，运行报告也需要在同一路径！！！
source  ~/app/bioinfo/dna/new.rc
/mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/bin/nextflow-1 run -bg /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/tdna.nf --fq_list $workdir/fq.list(rawdata) --ref ${ref.fa} --insert ${insert.fa} --outdir $workdir/output 
运行报告
运行nextflow给的outdir需要和git下来的dna文件夹在同一路径，运行报告也需要在同一路径！！！
source  ~/app/bioinfo/dna/new.rc
在运行环境创建data文件夹，data/project.info 和 data/info.log
这两个文件格式参考wgs结果02.reference同名文件，对应参考基因组信息查阅生信云即可。
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna_report.sh no ${tdna_结果路径}
bash /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/run_tdna_report.sh no /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/test/tang/raw_nf/output
```

## 原版
## 完成原始数据和基因组下载后，将TDNA的全序列获得，生成一个fasta文件，运行如下命令：

### TDNAscan
使用tdnascan运行，寻找单端read既比对到gdna，也比对到tdna上的结果
```
step1： tdnascan.py
        usage: tdnascan.py [-h] [--version] -1 FQ1 -2 FQ2 -g REFERENCE -t TDNA_SEQ -p
                   PROJECT [-n MINRD] [-a WINCLR] [-b WINDIR] [-@ THREAD]
```

使用aimhii运行，寻找read1和read2分别比对到tdna和gdna上的结果
```
step2： docker run -v "/mnt/ilustre:/mnt/ilustre" -w $PWD -i 10.2.4.236:5000/aimhii
        usage
            aimhii GENOME INSERT ADAPTER FASTQ1 FASTQ2 --outfile results.csv --plot readplot
```        


### 其他

ADAPTER文件如下：
```
        >TruSeq_Adapter_Index_1
        GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
```
