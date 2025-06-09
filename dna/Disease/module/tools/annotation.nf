tempdir = params.outdir
db = params.db

/*
 * Process ChrRecode: 根据染色体过滤vcf并且重命名染色体
 */
process chr_recode {
  publishDir "$tempdir/vcf", pattern: "raw.chr_recode.vcf.gz"
  tag "chr_recode"
  cpus 8
  memory "32G"
  input:
    path vcf
    path ped
    val rename

  output:
    path "raw.chr_recode.vcf.gz", emit: outvcf

  script:
  if(rename)
    """
    tabix -f ${vcf}
    cut -f 2 ${ped} > samples.txt
    bcftools view -S samples.txt --threads 8 -Oz -o chr.vcf.gz -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,sca1 ${vcf}
    # bcftools view --threads 8 -Oz -o chr.vcf.gz -t chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,sca1 ${vcf}
    bcftools annotate --threads 8 --rename-chrs ${projectDir}/bin/rename_chr.txt -o raw.chr_recode.vcf.gz -Oz chr.vcf.gz
    """
  else
    """
    tabix -f ${vcf}
    bcftools view --threads 8 -Oz -o raw.chr_recode.vcf.gz -t ${chrs} ${vcf}
    """
}

process sex_check {
  publishDir "$tempdir/annotation"
  tag "sex_check"
  cpus 8
  memory "32G"
  input:
    path vcf
    path ped

  output:
    path "sex_check.xls"

  script:
    """
    plink --threads 8 --vcf ${vcf} --split-x hg38 no-fail --make-bed --out testsex --double-id
    plink --threads 8 --check-sex --bfile testsex --fam ${ped} --out testsex
    cat testsex.sexcheck |tr -s ' ' '\\t' |sed 's/^\\t//g' > sex_check.xls
    """
}


/*
 * Process Coseparation: 共分离分析
 * TODO: 检查多等位新生突变的结果
 */
process coseparation_denovo {
  publishDir "$tempdir/vcf", pattern: "cog_denovo.vcf.gz"
  tag "coseparation"
  cpus 8
  memory "80G"
  input:
    path vcf
    path ped
    val cogped
  output:
    path "cog_denovo.vcf.gz", emit: outvcf
  script:
  if(cogped != "")
    """
    source ~/app/bioinfo/dna/miniconda3/etc/profile.d/conda.sh
    conda activate genmod
    gunzip -c ${vcf} | sed "s/Number='\\.'/Number=./g" | genmod models - --processes 8 --keyword ANN -t ped -f ${cogped} > cog.vcf
    bgzip -@ 8 cog.vcf
    tabix cog.vcf.gz
    triodenovo --ped ${cogped} --in_vcf cog.vcf.gz --out_vcf denovo.vcf
    bcftools query -f "%CHROM\\t%POS\\tDenovo_Mutation\\n" denovo.vcf > denovo.txt
    bgzip denovo.txt
    tabix -s 1 -b2 -e2 denovo.txt.gz
    bcftools annotate --threads 8 -a denovo.txt.gz -c CHROM,POS,Denovo_Mutation -h ${projectDir}/bin/denovo.header -Oz -o cog_denovo.vcf.gz cog.vcf.gz
    """
  else
    """
    gunzip -c ${vcf} | sed "s/Number='\\.'/Number=./g" | bgzip -c > cog.vcf.gz
    tabix cog.vcf.gz
    touch denovo.txt
    bgzip denovo.txt
    tabix -s 1 -b2 -e2 denovo.txt.gz
    bcftools annotate --threads 8 -a denovo.txt.gz -c CHROM,POS,Denovo_Mutation -h ${projectDir}/bin/denovo.header -Oz -o cog_denovo.vcf.gz cog.vcf.gz
    """
}


/*
 * Process BcftoolsNorm: 拆分多等位用于注释
 */
process bcftools_norm {
  publishDir "$tempdir/vcf", pattern: "normalized.vcf.gz"
  tag "bcftools_norm"
  cpus 8
  memory "32G"
  input:
    path vcf

  output:
    path "normalized.vcf.gz", emit: outvcf

  script:
    """
    bcftools norm --threads 8 -f ${db}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --multiallelics - --check-ref e --old-rec-tag Origin -Oz -o normalized.vcf.gz ${vcf}
    """
}



/*
 * Process FreqAnnot: 注释人群频率
 */
process freq_anno {
  publishDir "$tempdir/vcf"
  tag "freq_anno"
  cpus 8
  memory "64G"
  input:
    path vcf
    val type

  output:
    path "freq.vcf.gz", emit: outvcf

  script:
    """
    tabix ${vcf}
    bcftools annotate --threads 8 -a ${db}/gnomad.${type}s.hg38.txt.gz -C ${db}/gnomad.exomes.hg38.columns -h ${db}/gnomad.exomes.hg38.header -o gnomad_annotated.vcf.gz ${vcf}
    tabix gnomad_annotated.vcf.gz
    bcftools annotate --threads 8 -a ${db}/1000g.hg38.txt.gz -C ${db}/1000g.hg38.columns -h ${db}/1000g.hg38.header -o 1000g.vcf.gz gnomad_annotated.vcf.gz
    tabix 1000g.vcf.gz
    bcftools annotate --threads 8 -a ${db}/esp6500.hg38.txt.gz -C ${db}/esp6500.hg38.columns -h ${db}/esp6500.hg38.header -o esp6500.vcf.gz 1000g.vcf.gz
    tabix esp6500.vcf.gz
    bcftools annotate -a ${db}/dbsnp_156.vcf.gz -c ID -Oz -o freq.vcf.gz esp6500.vcf.gz
    rm gnomad_annotated.vcf.gz*
    rm 1000g.vcf.gz*
    rm esp6500.vcf.gz
    """
}


/*
 * Process NsfpAnnot: 注释位点功能
 */
process nsfp_anno {
  publishDir "$tempdir/vcf", pattern: "final*"
  tag "nsfp_anno"
  cpus 8
  memory "64G"
  input:
    path vcf

  output:
    path "final_region.vcf.gz", emit: outvcf

  script:
    """
    tabix ${vcf}
    bcftools annotate --threads 8 -a ${db}/dbscSNV1.1.hg38.txt.gz -C ${db}/dbscSNV1.1.hg38.columns -h ${db}/dbscSNV1.1.hg38.header -Oz -o sc.vcf.gz ${vcf}
    tabix sc.vcf.gz
    bcftools annotate --threads 8 -a ${db}/gerp++gt2.hg38.sort.txt.gz -C ${db}/gerp++gt2.hg38.sort.columns -h ${db}/gerp++gt2.hg38.sort.header -Oz -o gerp.vcf.gz sc.vcf.gz
    tabix gerp.vcf.gz
    bcftools annotate --threads 8 -a ${db}/gwas_catalog_v1.0.hg38.sorted.txt.gz -C ${db}/gwas_catalog_v1.0.hg38.sorted.txt.columns -h ${db}/gwas_catalog_v1.0.hg38.sorted.txt.header -Oz -o gwas.vcf.gz gerp.vcf.gz
    tabix gwas.vcf.gz
    bcftools annotate --threads 8 -a ${db}/dbNSFP4.4a.hg38.txt.gz -C ${db}/dbNSFP4.4a.columns -h ${db}/dbNSFP4.4a.hg38.header -Oz -o final.vcf.gz gwas.vcf.gz
    tabix final.vcf.gz
    bcftools annotate --threads 8 -a ${db}/targetScanS.hg38.sorted.bed.gz -C ${db}/targetScanS.hg38.sorted.columns -h ${db}/targetScanS.hg38.sorted.header -Oz -o final.tss.vcf.gz final.vcf.gz
    tabix final.tss.vcf.gz
    bcftools annotate --threads 8 -a ${db}/tfbsConsSites.hg38.sorted.bed.gz -C ${db}/tfbsConsSites.hg38.sorted.columns -h ${db}/tfbsConsSites.hg38.sorted.header -Oz -o final.tcs.vcf.gz final.tss.vcf.gz
    tabix final.tcs.vcf.gz
    bcftools annotate --threads 8 -a ${db}/hgmd.hg38.txt.gz -C ${db}/hgmd.hg38.txt.columns -h ${db}/hgmd.hg38.txt.header -Oz -o final.hgmd.vcf.gz final.tcs.vcf.gz
    tabix final.hgmd.vcf.gz
    bcftools annotate --threads 8 -a ${db}/repeatmask.hg38.sorted.txt.gz -C ${db}/repeatmask.hg38.sorted.columns -h ${db}/repeatmask.hg38.sorted.header -Oz -o final_region.vcf.gz final.hgmd.vcf.gz
    rm sc.vcf.gz*
    rm gerp.vcf.gz*
    rm gwas.vcf.gz*
    rm final.vcf.gz*
    rm final.tss.vcf.gz*
    rm final.tcs.vcf.gz*
    rm final.hgmd.vcf.gz*
    """
}


/*
 * Process AnnovarAnnot: annovar注释位点功能
 */
process annovar_anno {
  publishDir "$tempdir/vcf", pattern: "*.xls"
  tag "annovar_anno"
  cpus 8
  memory "64G"
  input:
    path vcf

  output:
    path "annovar.hg38_multianno.xls", emit: table

  script:
    """
    tabix ${vcf}
    bcftools query --allow-undef-tags -f "%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%ID\\t%INFO/AF_1KG_EAS\\t%INFO/AF_1KG\\t%INFO/AF_ESP6500\\t%INFO/AF_GNOMAD\\t%INFO/AF_GNOMAD_EAS\\t%INFO/dbscSNV_ADA_score\\t%INFO/dbscSNV_RF_score\\t%INFO/dbNSFP_Ensembl_transcriptid\\t%INFO/dbNSFP_Ensembl_proteinid\\t%INFO/dbNSFP_SIFT_score\\t%INFO/dbNSFP_SIFT_pred\\t%INFO/dbNSFP_Polyphen2_HDIV_score\\t%INFO/dbNSFP_Polyphen2_HDIV_pred\\t%INFO/dbNSFP_Polyphen2_HVAR_score\\t%INFO/dbNSFP_Polyphen2_HVAR_pred\\t%INFO/dbNSFP_LRT_score\\t%INFO/dbNSFP_LRT_pred\\t%INFO/dbNSFP_MutationTaster_score\\t%INFO/dbNSFP_MutationTaster_pred\\t%INFO/dbNSFP_MutationAssessor_score\\t%INFO/dbNSFP_MutationAssessor_pred\\t%INFO/dbNSFP_FATHMM_score\\t%INFO/dbNSFP_FATHMM_pred\\t%INFO/dbNSFP_PROVEAN_score\\t%INFO/dbNSFP_PROVEAN_pred\\t%INFO/dbNSFP_MetaSVM_score\\t%INFO/dbNSFP_MetaSVM_pred\\t%INFO/dbNSFP_CADD_raw\\t%INFO/dbNSFP_CADD_phred\\t%INFO/dbNSFP_phyloP100way_vertebrate\\t%INFO/dbNSFP_phyloP470way_mammalian\\t%INFO/dbNSFP_SiPhy_29way_logOdds\\t%INFO/dbNSFP_Interpro_domain\\t%INFO/GERPgt2\\t%INFO/TargetScanS\\t%INFO/TfbsConsSites\\t%INFO/HGMD_URL\\t%INFO/RepeatMask\\t%INFO/GWAS_disease_Pubmed_pValue\\t%INFO/Origin\\t%INFO/GeneticModels\\t%INFO/Denovo_Mutation[\\t%GT:%AD:%DP:%GQ]\\n" ${vcf} |perl -ne 'chomp;@a=split(/\\t/,\$_);if(length(\$a[3]) <= length(\$a[4])){print \$_,"\\n";next}else{print join("\\t",@a[0..1],\$a[2]+length(\$a[3])-length(\$a[4]),@a[3..\$#a]),"\\n"}'> annovar.input
    ${db}/table_annovar.pl annovar.input ${db} -buildver hg38 -out annovar -protocol ensGene,clinvar_20221231,cytoBand,cpgIslandExt,genomicSuperDups,mcap,revel,wgRna,intervar -operation g,f,r,r,r,f,f,r,f -nastring . --thread 16 --otherinfo
    bcftools query -l ${vcf} | sed 's/^/sample_/g' |cat ${projectDir}/bin/annovar.header - | tr '\\n' '\\t' |sed 's/\\t\$/\\n/g' >annovar.header
    cat annovar.header <(sed '1d' annovar.hg38_multianno.txt) > annovar.hg38_multianno.xls
    """
}


/*
 * Process FinalStat: 整理注释结果
 */
process final_stat {
  publishDir "$tempdir/annotation", pattern: "*.{xls,pdf,png}"
  tag "final_stat"
  cpus 8
  memory "64G"
  input:
    path multianno
    path ped
    val disease
    val gene
    val cogped
  output:
    path "filtered.Deleterious.xls", emit: del
    path "filtered.ExonicFunc.xls", emit: exon
    path "filtered.Freq.xls", emit: freq
    path "filtered.Func.xls", emit: func
    path "*_stat.xls", emit: filter_stat
    path "filtered.gene.list", emit: gene
    path "dgn.list", emit: dgnlist
    path "filtered.Freq.var.list", emit: variant
    path "*.enrichment.*", emit: enrichout
    val true, emit: finished
script:
    if (gene != '') {
        if (cogped != '') {
            """
            Rscript ${projectDir}/bin/exome_variant_filter.R -i ${multianno} -o filtered -a ${db} --ped ${ped} -d '${disease}' -g ${gene} --cogped ${cogped}
            cat filtered.gene.list ${gene} | sort | uniq > filtered.gene.list.tmp
            mv filtered.gene.list.tmp filtered.gene.list
            """
        } else {
            """
            Rscript ${projectDir}/bin/exome_variant_filter.R -i ${multianno} -o filtered -a ${db} --ped ${ped} -d '${disease}' -g ${gene}
            cat filtered.gene.list ${gene} | sort | uniq > filtered.gene.list.tmp
            mv filtered.gene.list.tmp filtered.gene.list
            """
        }
    } else {
        if (cogped != '') {
            """
            Rscript ${projectDir}/bin/exome_variant_filter.R -i ${multianno} -o filtered -a ${db} --ped ${ped} -d '${disease}' --cogped ${cogped}
            """
        } else {
            """
            Rscript ${projectDir}/bin/exome_variant_filter.R -i ${multianno} -o filtered -a ${db} --ped ${ped} -d '${disease}'
            """
        }
    }
}

/*
 * Process disgenet: disgenet 疾病关联基因
 */
process disgenet {
  publishDir "$tempdir/annotation", pattern: "*"
  tag "disgenet"
  cpus 2
  memory "4G"
  input:
    path disease
    path gene
    path variant
    val number

  output:
    path "dgn.*.result.xls", optional: true, emit: table
    path "dgn.result.png", optional: true, emit: png
    path "dgn.result.pdf", optional: true, emit: pdf
    val true, emit: finished

  script:
    if(number.size() > 0)
      """
      head -1 ${db}/disgenet/disgenet_gene_disease.txt > dgn.gene.result.xls
      awk -v FS="\\t" -v OFS="\\t" 'FNR==NR{a[\$1];next}(\$4 in a){print}' ${gene} ${db}/disgenet/disgenet_gene_disease.txt | awk -v FS="\\t" -v OFS="\\t" 'FNR==NR{a[\$1];next}(\$1 in a){print}' ${disease} - >> dgn.gene.result.xls
      head -1 ${db}/disgenet/disgenet_variant_disease.txt > dgn.variant.result.xls
      awk -v FS="\\t" -v OFS="\\t" 'FNR==NR{a[\$1];next}(\$3 in a){print}' ${variant} ${db}/disgenet/disgenet_variant_disease.txt | awk -v FS="\\t" -v OFS="\\t" 'FNR==NR{a[\$1];next}(\$1 in a){print}' ${disease} - >> dgn.variant.result.xls
      Rscript ${projectDir}/bin/plot_disgenet.R --genefile dgn.gene.result.xls --varfile dgn.variant.result.xls --outfile dgn.result
      """
}


/*
 * Process phenolyzer: phenolyzer 关联基因
 */
process phenolyzer {
  publishDir "$tempdir/annotation", pattern: "*"
  tag "phenolyzer"
  cpus 8
  memory "40G"
  input:
    path disease
    path gene
    val number

  output:
    path "phenolyzer.*.{xls,png,pdf}", optional: true, emit: table
    val true, emit: finished

  script:
    if(number.size() > 0)
      """
      cut -f 2 ${disease} | split -l 1 -d - disease_
      ls disease_* > disease.list
      xargs -P 8 -a disease.list -i disease_annotation.pl -f {} --gene ${gene} -p -ph -logistic -out phenolyzer.{} -addon_gg DB_MENTHA_GENE_GENE_INTERACTION -addon_gg_weight 0.05
      xargs -P 8 -a disease.list -i Rscript ${projectDir}/bin/plot_phenolyzer.R --infile phenolyzer.{}.final_gene_list --outfile phenolyzer.{} --dfile {}
      """
}

/*
 * Process genemania: genemania 蛋白互作分析
 */
process genemania {
  publishDir "$tempdir/annotation", pattern: "*"
  tag "genemania"
  cpus 8
  memory "80G"
  input:
    path gene

  output:
    path "genemania.*.xls", emit: table
    val true, emit: finished

  script:
    """
    java -Xmx8G -cp ${db}/genemania-cy3-3.5.3.jar org.genemania.plugin.apps.GeneSanitizer --data ${db}/gmdata-2021-04-29 --organism 9606 ${gene} | sed '1d' > gene_mania.map
    cut -f 2 gene_mania.map | tr -s "\\n" > genemania.list
    rm -fr queries
    mkdir queries
    split -d -a 3 -l 100 genemania.list genemania.split
    mv genemania.split* queries/
    cd queries
    for i in genemania.split*; do
        echo "H. Sapiens" > \${i}.query
        cat \${i} | tr "\\n" "\\t" >> \${i}.query
        echo -e "\\ncoexp\\tcoloc\\tpi\\tgi\\tpath\\n50" >> \${i}.query
    done
    ls genemania.split*.query | xargs -P 4 -i java -Xmx10G -cp ${db}/genemania-cy3-3.5.3.jar org.genemania.plugin.apps.QueryRunner --data ${db}/gmdata-2021-04-29 --out flat --threads 1 {}
    for i in genemania*.query-results.report.txt; do
        cat \${i} | perl ${projectDir}/bin/extract_mania.pl > \${i/query-results.report.txt/final.result}
    done
    csvtk -Tt concat genemania.*.final.result > ../genemania.result.xls
    cd ..
    Rscript ${projectDir}/bin/plot_genemania.R --infile genemania.result.xls --outfile genemania.stat.xls
    """
}

