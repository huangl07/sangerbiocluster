  ### 组装
  hifiasm -o sample.asm -t32 -z20 --dual-scaf --hg-size 600m -s 0.4 m84178_231212_141242_s3.hifi_reads.bc1010.fastq  2> log
  awk '/^S/{print ">"$2;print $3}' sample.asm.hic.gfa > sample.asm.hic.fa
  minimap2 -t 16 -ax map-hifi sample.asm.hap1.fa m84178_231212_141242_s3.hifi_reads.bc1010.fastq|samtools view -hF 256 -|samtools sort -@ 8 -m 1G -o aligned.bam -T tmp.ali
  samtools index aligned.bam

  #######################
  purge_haplotigs readhist -b aligned.bam -g sample.asm.hic.hap1.p_ctg.fasta -t 8
  purge_haplotigs contigcov  -i aligned.bam.gencov -low 35 -h 175 -m 95
  purge_haplotigs purge -g sample.bp.hap1.fa  -c coverage_stats.csv
  bwa index curated.fasta
  bwa mem -5SP -t 8 curated.fasta HIC.1.fastq.gz HIC.2.fastq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o HiC.curated.bam
  ######################
  split_fa sample.asm.hic.hap1.p_ctg.fasta > sample.asm.hic.hap1.p_ctg.fasta.split
  minimap2 -t 16 -xasm5 -DP sample.asm.hic.hap1.p_ctg.fasta.split sample.asm.hic.hap1.p_ctg.fasta.split| gzip -c - > sample.asm.hic.hap1.p_ctg.fasta.split.self.paf.gz
  samtools view -h aligned.bam|paftools.js sam2paf - > aligned.paf
  pbcstat aligned.paf
  calcuts PB.stat > cutoffs 2>calcults.log
  purge_dups -2 -T cutoffs_manual -c PB.base.cov sample.asm.hic.hap1.p_ctg.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
  get_seqs -e dups.bed sample.asm.hic.hap1.p_ctg.fasta
  bwa index purged.fa
  bwa mem -5SP -t 8 purged.fa HIC.1.fastq.gz HIC.2.fastq.gz|samblaster|samtools view - -@ 14 -S -h -b -F 3340 -o HiC.purged.bam
  #########################
  


  ~/app/bioinfo/dna/HapHiC-main/utils/filter_bam.py  --NM 3 --threads 8  HiC.bam 1 |samtools view - -b -@ 14 -o HiC.filtered.bam
  ~/app/bioinfo/dna/HapHiC-main/haphic pipeline curated.haplotigs.fasta HiC.filtered.bam 21 --RE "GATC"
 
