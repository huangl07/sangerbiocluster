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
    --compare   <file>  input compare list
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
compare_file = file(params.compare)
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
    script:
    """
		source ~/app/bioinfo/dna/new.rc
        bcftools filter --threads 8  -O z -i "F_Missing <=${params.miss} && MAF > ${params.maf} && FORMAT/DP < ${params.dep}" ${vcf}  > pop.filtered.vcf.gz
        tabix pop.filtered.vcf.gz
        bcftools filter --threads 8 -O z -i "F_Missing <=0" pop.filtered.vcf.gz > pop.nomissing.vcf.gz
        tabix pop.nomissing.vcf.gz
        bcftools view  -h  pop.filtered.vcf.gz |grep "contig=" |sed "s/##contig=<ID=//" |sed "s/,/\t/g" |sed "s/length=//g" |sed "s/>//g" |cut -f 1 > nomissing.list
        bcftools view  -h  pop.filtered.vcf.gz |grep "contig=" |sed "s/##contig=<ID=//" |sed "s/,/\t/g" |sed "s/length=//g" |sed "s/>//g" |cut -f 1 > chr.list
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
	source ~/app/bioinfo/dna/new.rc
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
        file compare from compare_file
   output:
        file "compare.list" into compare,compare1
   script:
   """
    perl ${baseDir}/bin/xpehh.pl -d ./ -o compare.list -c ${compare}
   """
}
compare.splitCsv(header:false,sep:'\t').map{row-> tuple(file(row[0]),file(row[1]),row[2],row[3])}.set{xpehh}
process selscan{
    publishDir "${params.outdir}/12.selscan", pattern:"*"
   queue "SANGERDEV"
   cpus 8
   executor "slurm"
   memory "30G"
   input:
      tuple file(g1),file(g2),compare,chr  from xpehh
   output:
       file "${compare}-${chr}.xpehh.out" into  xpehh_result
   script:
   """
   source ~/app/bioinfo/dna/new.rc
    plink --vcf ${g1} --recode --allow-extra-chr --double-id --chr-set 55
    Rscript ${baseDir}/bin/map-xpehh.R --map plink.map --out total.map
    selscan --xpehh --vcf ${g1} --vcf-ref ${g2}  --map total.map --threads 8 --out ${compare}-${chr} --unphased --cutoff 0.05 --trunc-ok
   """
}


process rearrange_xpehh_result{
    publishDir "${params.outdir}/13.xpehh", pattern:"*xpehh_abstract"
    queue "SANGERDEV"
    cpus 1
    executor "slurm"
    memory "30G"
    cache 'lenient'
    input:
        file xpehh from xpehh_result.collect()
    output:
        path "total.xpehh.windows.result" into xpehh_ch
        path "total.xpehh.norm.result" into norm_ch
    script:
    """
	source ~/app/bioinfo/dna/new.rc
    norm --xpehh  --files *.xpehh.out --winsize 10000 --qbins 1 --min-snps 1 --bp-win
    for i in *.xpehh.out.norm.*.windows; do
        chr=\${i/*-/}
        chr=\${chr/.*/}
        name=\${i/-*/}
        sed -i -e 's/^/'\${name}\\\\t\${chr}'\t/g' \$i
    done
    for i in *.xpehh.out.norm; do
        chr=\${i/*-/}
        chr=\${chr/.*/}
        name=\${i/-*/}
        sed  -i -e 's/^\\./'\${name}\\\\t\${chr}'/g' \$i
    done
    cat *.norm|grep -v "id" > total.xpehh.norm.result
    cat *.windows|grep -v "id" > total.xpehh.windows.result
    sed -i '1i Compare\\tChr\\tPos\\tGPos\\p1\\tihh1\\tp2\\tihh2\\txpehh\\tnormxpehh\\tcrit'   total.xpehh.norm.result  
    sed -i '1i Compare\\tChr\\tStart\\tEND\\tNSites\\tfrac1\\tfrac2\\tpercentile1\\tpercentile2\\tmax\\tmin'   total.xpehh.windows.result  
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
        path xpehh from xpehh_ch.flatten()
        path norm_ch from norm_ch.flatten()
    output:
        file "*"
        file "*.select" into select_ch4
    script:
    """
		source ~/app/bioinfo/dna/new.rc
        Rscript ${baseDir}/bin/xpehh-window.R --xpehh ${xpehh} --threshold 0.05
    """
}
















