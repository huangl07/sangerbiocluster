less -S germline.sniffles.anno.vcf.gz |grep -v "##"|less -S |perl -ne 'chomp;@a=split;if(/\#CHROM/){$id1=$a[9];$id2=$a[10];next;}next if(!/INS/);$b=(split(":",$a[9]))[0];$c=(split(":",$a[10]))[0];next if ($a[4] =~ /[^ATGCN]/);print join("\n",">$a[0]"."-".$a[1]."-".$id1."-".$b."-".$id2."-".$c,$a[4]),"\n"'|less -S > insert.fa
makeblastdb -in insert.fa -dbtype nucl
 blastn -query tdna.fa -db insert.fa -num_threads 8 -out blast.out  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"
  less -S blast.out |perl -ne '@a=split;@b=split(/\-/,$a[1]);if(!defined $nohead){print "#seqid\tstart\tend\tchromosome\tinsert_pos\t$b[2]\t$b[4]\tinsert_start\tinsert_end\tinsert_lenth\ttdna_lenth\n";$nohead=1}print join("\t",$a[0],$a[6],$a[7],$b[0],$b[1],$b[3],$b[5],$a[8],$a[9],$a[13],$a[14]),"\n";' > result.table.xls


#seqid: 带转入序列id
#start: 转入基因组的序列起点
#end：转入基因组的序列终点
#chromosome: 转入基因组的染色体编号
#insert_pos: 转入基因组的位置
#A12K：A12K在该位置的基因型，0/0为和基因组一致，1/1为纯和，0/1为杂合
#D31：与A12K一致
#insert_start：转基因序列比对insert序列的起点
#insert_end: 转基因序列比对insert序列的终点
#insert_lenth：insert序列的总长度
#tdna_lenth：tdna序列的长度
#
#方法：使用minimap2将ont测序数据比对参考基因组后，使用sniffiles进行svcalling，calling后提取所有的insert变异序列，使用blast软件将insert序列与TDNA序列进行比对，比对到的就是TDNA插入的位置
#
