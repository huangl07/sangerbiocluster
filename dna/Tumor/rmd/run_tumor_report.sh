#!bin/bash

task_id=$1
workdir=`pwd`
data_realease_dir="$workdir/data_release"
result_dir=$2
workflow_results=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/`
echo $workflow_results
resultdir=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/tmp`
filedir=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id/output/published/data`
#WES报告
cd $workdir

if [ ! -d $data_realease_dir ];then
    mkdir -p $data_realease_dir
fi

## 报告模板导入
if [ ! -d ./report ];then
    mkdir -p ./report/Result
    mkdir -p ./report/file
    mkdir -p ./report/info
fi


##sample_tumor pair处理
pair_info=$3
awk '{print $1"_"$2}' $pair_info > ./report/info/pair.info

cp -r /mnt/lustre/users/sanger-dev/sg-users/yexiaoya/human_exon/hrda/rmd/ $workdir/report

#WES报告
info_dir="$workdir/report/info"

cp $resultdir/02.reference/project.info $info_dir/project.info
cp $resultdir/02.reference/info.log $info_dir/info.log
cp $resultdir/01.fastq_qc/qc.stat $info_dir
cp $resultdir/04.snpIndel/snp/snp.stat.all $info_dir
cp $resultdir/04.snpIndel/snp/snp.stat $info_dir
cp $resultdir/04.snpIndel/snp/snp_anno.stat $info_dir
cp $resultdir/04.snpIndel/indel/indel.stat $info_dir
cp $resultdir/01.fastq_qc/json/* $info_dir
cp $resultdir/03.mappingStat/all.summary.stats $info_dir
cp $resultdir/03.mappingStat/*.png $info_dir
cp $resultdir/08.circos/circos.png $info_dir
#等上线以后更改
cp $resultdir/01.fastq_qc/bammatcher/bam_matcher.report $info_dir/bam_matcher.report
patient_num=$(cut -f1 --complement $workdir/data_info/info.txt | uniq  | grep -v "-" | wc -l)
sample_num=$(cut -f2 --complement $workdir/data_info/info.txt | uniq  | grep -v "-" | wc -l)
control_num=$(cut -f3 --complement $workdir/data_info/info.txt | uniq  | grep -v "-" | wc -l)
echo -e "${patient_num}\t${sample_num}\t${control_num}" > $info_dir/info.stat

awk '{print $1,$2,$3,$4,$5}' $resultdir/04.snpIndel/indel/indel.stat | grep -v "pop" |sed 's/ /\t/g'|sed 1d |sed '1i\Sample ID\tInsertion Number\tDeletion Number\tHeterozygosity Number\tHomozygosity Number'>$info_dir/indel_stat.xls
cut -f1,16,17,20,27,29,30 $resultdir/04.snpIndel/indel/indel_anno.stat |sed 1d |sed '1i\Sample ID\tExon Loss Variant\tFrameshift Variant\tIntragenic Variant\tStart Lost\tStop Gained\tStop Lost'>$info_dir/indel_anno.xls
cut -f1,6,7,8,9 $resultdir/04.snpIndel/indel/indel_anno.stat |sed 's/sampleID/sample ID/g' >$info_dir/indel_anno_stat.xls
if [ -d $resultdir/05.sv ];then
    cp $resultdir/05.sv/stat/stat.txt $info_dir/sv_stat.xls
fi
if [ -d $resultdir/06.cnv ];then
    cp $resultdir/06.cnv/*.stats.xls $info_dir/cnv_stat.xls
fi
if [ -d $resultdir/07.ssr ];then
    cat $resultdir/07.ssr/stat/table_stat.out |sed 's/Sample/Sample ID/g' |sed 's/SSR_Number/SSR Number/g' >$info_dir/ssr_stat.xls
else
    cat ~/app/database/wgs_v4_bed/pdf_report_base/wgs/data/ssr_stat.xls |sed 's/Sample/Sample ID/g' |sed 's/SSR_Number/SSR Number/g' > $info_dir/ssr_stat.xls
fi
cat $resultdir/09.anno/pop.stat.xls  |sed 's/#type/type/g' >$info_dir/gene_anno.xls 

#sampleID:
grep -v "Duplication" $resultdir/06.cnv_somatic/all.cnvkit.stats.xls |cut -f1 >$info_dir/somatic_stat_1
#snv:
for j in `ls $resultdir/04.snpIndel_somatic/finalvcf/*somatic_snp_filtered_anno.vcf.gz `;do wc -l $j|cut -d " " -f1;done >$info_dir/somatic_stat_2
# wc -l $resultdir/04.snpIndel_somatic/finalvcf/*somatic_snp_filtered_anno.vcf.gz |grep -v "total" |cut -d ' ' -f4 >$info_dir/somatic_stat_2
#indel:
# wc -l $resultdir/04.snpIndel_somatic/finalvcf/*somatic_indel_filtered_anno.vcf.gz |grep -v "total" |cut -d ' ' -f4 >$info_dir/somatic_stat_3
for j in `ls $resultdir/04.snpIndel_somatic/finalvcf/*somatic_indel_filtered_anno.vcf.gz `;do wc -l $j|cut -d " " -f1;done >$info_dir/somatic_stat_3
#sv:
for j in `ls $resultdir/05.sv_somatic/stat/*/*stat.xls `;do tail -1 $j |awk '{print $2+$3+$4+$5+$6}'; done >$info_dir/somatic_stat_4
#cnv:
grep -v "Duplication" $resultdir/06.cnv_somatic/all.cnvkit.stats.xls|awk '{print $3+$4}' >$info_dir/somatic_stat_5
#总表
paste $info_dir/somatic_stat_1 $info_dir/somatic_stat_2 $info_dir/somatic_stat_3 $info_dir/somatic_stat_4 $info_dir/somatic_stat_5 >$info_dir/somatic_stat.xls1
sed '1i\Tumor Normal\tSNV Number\tIndel Number\tSV Number\tCNV Number' $info_dir/somatic_stat.xls1 >$info_dir/somatic_stat.xls

#somatic-snp-anno
cat $resultdir/04.snpIndel_somatic/snp/*_snp_anno.stat |grep -v "3_prime_UTR_variant" |awk '{if($5!="0" || $6!="0" ||$7!="0" || $8!="0" ) print $0}' |cut -f10,12,17,18,19,21 >$info_dir/somatic.snp.anno.1
paste $info_dir/somatic_stat_1 $info_dir/somatic.snp.anno.1 >$info_dir/somatic.snp.anno.2
sed '1i\Sample ID\tIntergenic Region\tMissense Variant\tStart Lost\tStop Gained\tStop Lost\tSynonymous Variant' $info_dir/somatic.snp.anno.2 >$info_dir/somatic.snp.anno.xls

#somatic-indel-anno
cat $resultdir/04.snpIndel_somatic/indel/*_indel_anno.stat |grep -v "3_prime_UTR_variant" |awk '{if($9!="0" || $6!="0" ||$7!="0" || $8!="0" ) print $0}' |cut -f16,17,20,27,29,30 >$info_dir/somatic.indel.anno.1
paste $info_dir/somatic_stat_1 $info_dir/somatic.indel.anno.1 >$info_dir/somatic.indel.anno.2
sed '1i\Sample ID\tExon Loss Variant\tFrameshift Variant\tIntragenic Variant\tStart Lost\tStop Gained\tStop Lost'  $info_dir/somatic.indel.anno.2 >$info_dir/somatic.indel.anno.xls


#tumor报告
cd $data_realease_dir

mkdir 00.basic_anlysis    01.purity_ploidy    02.predisposing_gene    03.oncogene    04.clone_anlysis   05.CNV-LOH    06.mutation_spectrum    07.hypermutation_classification    08.significantly_mutations    09.hla_subtype    10.netMHC    11.drug_target

# 00.basic_anlysis结果整理
cp $result_dir/00.basic_anlysis/pop.snpindel.final.vcf ./00.basic_anlysis/pop.snpindel.final.vcf
cp $result_dir/00.basic_anlysis/tcga_mutect.result ./00.basic_anlysis/TCGA_mutect.vcf
mkdir -p 00.basic_anlysis/somatic_vcf
cp $result_dir/00.somatic_vcf/*.vcf ./00.basic_anlysis/somatic_vcf/

# 01.cnl结果整理
touch ./01.purity_ploidy/purity_ploidy.stat
echo "Tumor_control	Cellularity	Ploidy.estimate	Ploidy.mean.cn" >> ./01.purity_ploidy/purity_ploidy.stat
for i in $(cut -f1 $workdir/report/info/pair.info);do
  cp $result_dir/01.cnl/sequenza_plot/$i/${i}_CP_contours.pdf ./01.purity_ploidy/${i}_CP_contours.pdf
  cp $result_dir/01.cnl/sequenza_plot/$i/${i}_CP_contours.png ./01.purity_ploidy/${i}_CP_contours.png
  echo -e "${i}\t$(awk 'NR==2' $result_dir/01.cnl/sequenza_plot/$i/${i}_confints_CP.txt)" >> ./01.purity_ploidy/purity_ploidy.stat
done
cp ./01.purity_ploidy/purity_ploidy.stat $info_dir/sequenza.stat

# 02.pre_gene结果整理
touch ./02.predisposing_gene/pre_gene.stat
echo "Tumor_control	Gene_numbers" >> ./02.predisposing_gene/pre_gene.stat
for i in $(cut -f1 $workdir/report/info/pair.info);do
  name=``${i}``
  a=$(cut -f22 --complement $result_dir/02.pre_gene/CGC_anno/${i}_CGC_anno.txt | uniq  | grep -v "-" | wc -l)
  echo -e "${name}\t${a}" >> ./02.predisposing_gene/pre_gene.stat
  cp $result_dir/02.pre_gene/CGC_anno/${i}_CGC_anno.txt  ./02.predisposing_gene/${i}.merge_anno.txt 
  cp $result_dir/02.pre_gene/CGC_anno/${i}_CGC_anno.xls  ./02.predisposing_gene/${i}.merge_anno.xls 
done
cp ./02.predisposing_gene/pre_gene.stat $info_dir/pre_gene.stat

# 03.oncogene结果整理
touch ./03.oncogene/oncogene.stat
echo "Sample_Control	Genes_all	IntOGen	CGC	Com299	Com435	B125" >> ./03.oncogene/oncogene.stat
cp $result_dir/03.oncogene/* ./03.oncogene/
for i in $(cut -f1 $workdir/report/info/pair.info);do
  name=``${i}``
  a=$(cut -f1 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq  | grep -v "-" | wc -l)
  b=$(cut -f1,7 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq | grep -v "-" | wc -l)
  c=$(cut -f1,11 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq | grep -v "-" | wc -l)
  d=$(cut -f1,12 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq | grep -v "-" | wc -l)
  e=$(cut -f1,13 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq | grep -v "-" | wc -l)
  f=$(cut -f1,16 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq | grep -v "-" | wc -l)
  echo -e "${name}\t${a}\t${b}\t${c}\t${d}\t${e}\t${f}" >> ./03.oncogene/oncogene.stat >> ./03.oncogene/oncogene.stat
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}_plot_oncodriven.pdf ./03.oncogene/${i}_plot_oncodriven.pdf
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}_plot_oncodriven.png ./03.oncogene/${i}_plot_oncodriven.png
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}_plot_oncodriven.png $info_dir/${i}_plot_oncodriven.png
  cp $result_dir/03.oncogene/${i}_oncogene.txt ./03.oncogene/${i}_oncogene.txt 
  cp $result_dir/03.oncogene/${i}_oncogene.xls ./03.oncogene/${i}_oncogene.xls
done
cp ./03.oncogene/oncogene.stat $info_dir/oncogene.stat


# 04.clone_anlysis结果整理
for i in $(cut -f1 $workdir/report/info/pair.info);do
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}_oncogene.xls  ./04.clone_anlysis/${i}_oncogene.xls
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}.pdf  ./04.clone_anlysis/${i}.pdf
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}.png  ./04.clone_anlysis/${i}.png
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}.png  $info_dir/${i}_clone.png
  cp $result_dir/04.clone_anlysis/clone_analysis/${i}cluster_data.txt ./04.clone_anlysis/${i}_cluster_data.txt
done

# 05.CNV-LOH结果整理
for i in $(cut -f1 $workdir/report/info/pair.info);do
  cp $result_dir/01.cnl/sequenza_plot/$i/$i/ascat_hg38_cov_hist.pdf ./05.CNV-LOH/${i}_CNV_LOH.pdf
  cp $result_dir/01.cnl/sequenza_plot/$i/$i/ascat_hg38_cov_hist.pdf $info_dir/${i}_CNV_LOH.pdf
  convert $result_dir/01.cnl/sequenza_plot/$i/$i/ascat_hg38_cov_hist.pdf $result_dir/01.cnl/sequenza_plot/$i/$i/ascat_hg38_cov_hist.png
  cp $result_dir/01.cnl/sequenza_plot/$i/$i/ascat_hg38_cov_hist.png ./05.CNV-LOH/${i}_CNV_LOH.png
  cp $result_dir/01.cnl/sequenza_plot/$i/$i/ascat_hg38_cov_hist.png $info_dir/${i}_CNV_LOH.png
done

# 06.hf_mutation_gene_anlysis结果整理
mkdir -p 06.mutation_spectrum/mutation_spectrum_sample_control
mkdir -p 06.mutation_spectrum/mutation_spectrum_merge
cp $result_dir/05.hf_mutation_gene_anlysis/mutation_spectrum/* ./06.mutation_spectrum/mutation_spectrum_sample_control/
cp $result_dir/05.hf_mutation_gene_anlysis/mutation_spectrum_merge/* ./06.mutation_spectrum/mutation_spectrum_merge/
cp $result_dir/05.hf_mutation_gene_anlysis/mutation_spectrum_merge/all_mutation_pattern.png $info_dir/mutation_pattern_all.png
cp $result_dir/05.hf_mutation_gene_anlysis/mutation_spectrum_merge/all_mutation_spectrum.png $info_dir/mutation_spectrum_all.png

cp $result_dir/05.hf_mutation_gene_anlysis/mutation_spectrum/* $info_dir/

# 07.hypermutation_classification结果整理
mkdir -p ./07.hypermutation_classification/result/
cp $result_dir/06.hypermutation_classification/msi_result_all.xls ./07.hypermutation_classification/
cp $result_dir/06.hypermutation_classification/msi_result_all.txt ./07.hypermutation_classification/
cp $result_dir/06.hypermutation_classification/result/*.prefix_somatic ./07.hypermutation_classification/result/
cp $result_dir/06.hypermutation_classification/result/*.prefix_germline ./07.hypermutation_classification/result/
cp $result_dir/06.hypermutation_classification/result/*.prefix_all ./07.hypermutation_classification/result/
cp $result_dir/06.hypermutation_classification/result/*.prefix_unstable ./07.hypermutation_classification/result/
cp $result_dir/06.hypermutation_classification/msi_result_all.txt $info_dir/msi.stat

# 08.significantly_mutations结果整理
cp $result_dir/07.significantly_mutated/mutation_relation/panorama_of_genomic_mutations.png ./08.significantly_mutations/panorama_of_genomic_mutations.png
cp $result_dir/07.significantly_mutated/mutation_relation/panorama_of_genomic_mutations.png $info_dir/panorama_of_genomic_mutations.png
cp $result_dir/07.significantly_mutated/mutation_relation/panorama_of_genomic_mutations.pdf ./08.significantly_mutations/panorama_of_genomic_mutations.pdf
cp $result_dir/07.significantly_mutated/music2/smgs ./08.significantly_mutations/smgs.txt
cp $result_dir/07.significantly_mutated/music2/smgs $info_dir/smgs.txt
cp -r $result_dir/07.significantly_mutated/enrich/GO_result ./08.significantly_mutations/GO_result
cp -r $result_dir/07.significantly_mutated/enrich/KEGG_result ./08.significantly_mutations/KEGG_result
if [ -f $result_dir/07.significantly_mutated/enrich/KEGG_result/all_KEGGenrichment.png ]; then 
  cp $result_dir/07.significantly_mutated/enrich/KEGG_result/all_KEGGenrichment.png $info_dir/;
fi
if [ -f $result_dir/07.significantly_mutated/enrich/GO_result/all_GOenrichment.png ]; then 
  cp $result_dir/07.significantly_mutated/enrich/GO_result/all_GOenrichment.png $info_dir/;
fi
if [ ! -f $result_dir/07.significantly_mutated/mutation_relation/error.log ]; then 
  cp $result_dir/07.significantly_mutated/mutation_relation/mutation_relation.png ./08.significantly_mutations/;
  cp $result_dir/07.significantly_mutated/mutation_relation/mutation_relation.pdf ./08.significantly_mutations/;
  cp $result_dir/07.significantly_mutated/mutation_relation/mutation_relation.png $info_dir/mutation_relation.png
fi

# 09.hla_subtype 结果整理
cp $result_dir/08.hla_subtype/result/* ./09.hla_subtype/
cp $result_dir/08.hla_subtype/sample/* ./09.hla_subtype/
cut -f1,2,3,4,5 $result_dir/08.hla_subtype/result/MHC_result_all.txt > $info_dir/hla.stat

# 10.netMHC 结果整理
mkdir -p ./10.netMHC/MHC_I/
mkdir -p ./10.netMHC/MHC_II/
mkdir -p ./10.netMHC/MHC_II/detail/
mkdir -p ./10.netMHC/MHC_I/detail/
mkdir -p ./10.netMHC/MHC_I/final_stat/
mkdir -p ./10.netMHC/MHC_II/final_stat/
cp $result_dir/09.netMHC/MHC_I/final_xlsx/* ./10.netMHC/MHC_I/detail/
cp $result_dir/09.netMHC/MHC_I/final_stat_txt/* ./10.netMHC/MHC_I/final_stat/
cp $result_dir/09.netMHC/MHC_II/final_xlsx/* ./10.netMHC/MHC_II/detail/
cp $result_dir/09.netMHC/MHC_II/final_stat_txt/* ./10.netMHC/MHC_II/final_stat/
cp -r $result_dir/09.netMHC/MHC_I/final_stat_txt/* $info_dir/

# 11.drug_target 结果整理
cp $result_dir/10.drug_target/result/* ./11.drug_target/
touch ./11.drug_target/drug.stat
echo "Sample_Control	Gene_numbers	Drug_numbers" >> ./11.drug_target/drug.stat
for i in $(cut -f1 $workdir/report/info/pair.info);do
  name=``${i}``
  a=$(cut -f7 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq  | grep -v "-" | wc -l)
  b=$(cut -f10 --complement $result_dir/03.oncogene/${i}_oncogene.txt | uniq  | grep -v "-" | wc -l)
  echo -e "${name}\t${a}\t${b}" >> ./11.drug_target/drug.stat
done
cp ./11.drug_target/drug.stat $info_dir/drug.stat

### 把data_release的结果放到report/Result
cp $workdir/data_info/info.txt $workdir/report/info/info.txt
cd $workdir
cp -r $workdir/data_release $workdir/report/Result
cp -r $workdir/report/info $workdir/report/rmd/
cp $3 $workdir/report/rmd/info/tumor_pair_info.txt
cp $workdir/data_info/info.txt $workdir/report/rmd/info/info.txt
cp $workdir/report/info/info.stat $workdir/report/rmd/info/info.stat


if test -d "$resultdir/05.sv"; then
    svy="Yes"
else
    svy="No"
fi
echo ${svy}
if test -d "$resultdir/06.cnv"; then
    cnv="Yes"
else
    cnv="No"
fi
echo ${cnv}
Rscript $workdir/report/rmd/rmarkdown.R --rmd $workdir/report/rmd/report.rmd --run_sv ${svy} --run_cnv ${cnv} --format html --outfile $workdir/report/tumor_all.html
Rscript $workdir/report/rmd/rmarkdown.R --format pdf --outfile $workdir/report/tumor_all.pdf --rmd $workdir/report/rmd/report.rmd --run_sv ${svy} --run_cnv ${cnv}  