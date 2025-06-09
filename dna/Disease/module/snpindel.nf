include { chr_recode; coseparation_denovo; bcftools_norm; freq_anno; nsfp_anno; annovar_anno; disgenet; phenolyzer; sex_check; genemania; final_stat } from './tools/annotation'


workflow annot_snpindel {
  take:
    vcf
    rename
    type
    ped
    cog
    disease
    cogped
    gene
    noasso

  main:
    chr_recode(file(vcf), file(ped), rename)
    sex_check(chr_recode.out.outvcf, file(ped))
    if (cog) {
      coseparation_denovo(chr_recode.out.outvcf, file(ped), cogped)
      bcftools_norm(coseparation_denovo.out.outvcf)
    } else {
      bcftools_norm(chr_recode.out.outvcf)
    }
    freq_anno(bcftools_norm.out.outvcf, type)
    nsfp_anno(freq_anno.out.outvcf)
    annovar_anno(nsfp_anno.out.outvcf)
    final_stat(annovar_anno.out.table, file(ped), disease, gene, cogped)
    receiver = final_stat.out.finished.collect()
    if(!noasso){
      genemania(final_stat.out.gene)
      disgenet(final_stat.out.dgnlist, final_stat.out.gene, final_stat.out.variant, final_stat.out.dgnlist.splitText( by: 1 ).collect())
      phenolyzer(final_stat.out.dgnlist, final_stat.out.gene, final_stat.out.dgnlist.splitText( by: 1 ).collect())
      receiver = final_stat.out.finished.concat(genemania.out.finished, disgenet.out.finished, phenolyzer.out.finished).collect()
    }
    finished = receiver.collect()

  emit:
    finished = finished
}

