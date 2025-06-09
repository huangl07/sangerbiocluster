#!/usr/bin/env nextflow
params.outdir = "demo"
params.help = false
params.method="GATK"
params.miss=0.3
params.maf=0.05
params.win=1000000
params.dep=7
def helpMessage() {

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --fqlist fqlist --outdir '/project/'

    --vcf       <file>  input vcf file
    --outdir    <dir>   output dir
    --group     <file>  input group file
    --miss  <num>   miss
    --maf   <num>   maf
    --win   <num>   window size
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

vcf_file = files(params.vcf)
group_file = files(params.group)

process vcffilter{
    publishDir "${params.outdir}/01.filter", pattern:"*"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'

    input:
        file vcf from vcf_file
    output:
        file "pop.filtered.vcf.gz" into filter_vcf1,filter_vcf2,filter_vcf3,filter_vcf4,filter_vcf5,filter_vcf6
        file "pop.filtered.vcf.gz.tbi" into vcf_index
        file "pop.nomissing.vcf.gz" into nomissing
        file "nomissing.list" into nomissing_chr
        file "chr.list" into chr_list,chr_list1
        file "chr_length.list" into chr_length
    script:
    """
        bcftools filter --threads 8  -O z -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf.gz
        tabix pop.filtered.vcf.gz
        bcftools filter --threads 8 -O z -i "F_Missing <=0" pop.filtered.vcf.gz > pop.nomissing.vcf.gz
        tabix pop.nomissing.vcf.gz
        bcftools view  -h  pop.filtered.vcf.gz |grep "contig=" |sed "s/##contig=<ID=//" |sed "s/,/\t/g" |sed "s/length=//g" |sed "s/>//g" |cut -f1> nomissing.list
        bcftools view  -h  pop.filtered.vcf.gz |grep "contig=" |sed "s/##contig=<ID=//" |sed "s/,/\t/g" |sed "s/length=//g" |sed "s/>//g" |awk '{print \$1,int(\$2/10000)+1}' |sed 's/ /\t/'  > chr_length.list 
        cut -f1 chr_length.list >chr.list
    """
}


process group_preparie{
    publishDir "${params.outdir}/02.group", pattern:"*"
    queue "SANGERDEV"
    cpus 8
    executor "slurm"
    memory "30G"
        cache 'lenient'

    input:
        file group from group_file
    output:
        file "split.list" into grouplist1,grouplist2,tmp
        file "xpclr.list" into xpclr_group,xpehh_group
        file "treemix.list" into treemix_group,treemix1
    script:
    """
        perl ${baseDir}/bin/group.pl -i ${group} -o ./
        perl ${baseDir}/bin/remakegrolist.pl -in ${group} -out  treemix.list
     """
}
grouplist1.splitCsv(header:false,sep:'\t').map{row-> tuple(row[0],file(row[1]))}.combine(filter_vcf1).set{groups}
xpclr_group.splitCsv(header:false,sep:'\t').combine(filter_vcf3).set{groups3}

//####################################################################################

grouplist2.splitCsv(header:false,sep:'\t').map{row-> tuple(row[0],file(row[1]))}.combine(nomissing).set{group4}

nomissing_chr.splitCsv().combine(group4).set{chr_vcf_input}

process vcf_chr{
   publishDir "${params.outdir}/10.vcfsplit", pattern:"*"
   queue "SANGERDEV"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
        tuple chr,gid,gtxt,vcf from chr_vcf_input
      // each chr from chrs
   output:
       file "${chr}-${gid}.vcf.gz" into  vcf_chr
   script:
   """
   bcftools view -r ${chr} -S ${gtxt} -m2 -M 2  --threads 8 $vcf > ${chr}-${gid}.vcf.gz
   """
}
process compare{
   publishDir "${params.outdir}/11.xpehhcompare", pattern:"*"
   queue "SANGERDEV"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
        file group from group_file
        file vcf from vcf_chr.collect()
        file chr_length from chr_length
   output:
        file "compare.list" into compare,compare1
   script:
   """
    readlink -f *.vcf.gz | awk -F/ '{
        path = \$0;
        split(path, path_parts, "/");
        file_name = path_parts[length(path_parts)];
        split(file_name, name_parts, "-");
        chromosome = name_parts[1];
        sample_name1 = name_parts[2];
        split(sample_name1, name_parts1, ".vcf.");
        sample_name = name_parts1[1];
        print path " " "hhhhhhhh" " " sample_name " " chromosome
    }' |sed 's/ /\t/g' >compare.list1
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[\$1]=\$0;next} (\$4 in a){print \$0, a[\$4]}' ${chr_length}  compare.list1 |cut -f 1,2,3,4,6  >compare.list
   """

}
compare.splitCsv(header:false,sep:'\t').map{row-> tuple(file(row[0]),file(row[1]),row[2],row[3],row[4])}.set{xpehh}
process selscan{
    publishDir "${params.outdir}/12.selscan", pattern:"*"
   queue "SANGERDEV"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
      tuple file(g1),file(g2),compare,chr,grid  from xpehh
   output:
       file "SweeD_Report.${compare}-${chr}.10kb" into  clr_result, optional
   script:
   """
   source ~/app/bioinfo/dna/new.rc
    SweeD -name ${compare}-${chr}.10kb -input ${g1} -grid ${grid}  -minsnps 200 || echo "no clr result"

   """
}


process rearrange_clr_result{
    publishDir "${params.outdir}/13.ihs", pattern:"*"
    queue "SANGERDEV"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        file xpehh from clr_result.collect()
    output:
        path "total.clr.result" into clr_ch
    script:
    """
        source ~/app/bioinfo/dna/new.rc
        for file in *.10kb; do
            if [ -s "\$file" ]; then
                sed -i 1d \$file
                sed -i 1d \$file
                sed -i 1d \$file
                base_name=\$(basename "\$file")
                chr=\${base_name#*-}    # 去除开头直到 "-" 的部分
                chr=\${chr%%.*}           # 去除结尾直到 "." 的部分
                name=\${base_name#*SweeD_Report.}
                name=\${name%%-*}
                sed  -i -e 's/^/'\${name}\\\\t\${chr}'\\t/g' \$file
            fi
        done
        cat *.10kb | awk '{print \$1,\$2,int(\$3),int(\$3),\$4,\$5,int(\$6),int(\$7)}' |sed 's/ /\\t/g'| sed '1i Compare\\tChromosome\\tPosition1\\tPosition2\\tLikelihood\\tAlpha\\tStartPos\\tEndPos'> total.clr.result
    """
}


process draw_manhattan_xpehh{
    publishDir "${params.outdir}/14.xpehhResult", pattern:"*"
    queue "SANGERDEV"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        path clr from clr_ch.flatten()
    output:
        file "*"
        file "*.select" into select_ch4
    script:
    """
        source ~/app/bioinfo/dna/new.rc
        Rscript ${baseDir}/bin/clr-window.R --xpehh ${clr} --threshold 0.05
    """
}
















