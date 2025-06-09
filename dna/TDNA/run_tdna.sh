source  /mnt/lustre/users/sanger-dev/app/bioinfo/dna/new.rc

#task_id=$1
workdir=`pwd`
#workflow_results=`readlink -f /mnt/ilustre/isanger_workspaceWgsV4/*/WgsV4_$task_id`
workflow_results=$1
DIRname=$(dirname $0)
if [ ! -f "$DIRname/dna/TDNA/tdna_wgs.nf" ]; then
  echo "No GIT !!! Using default tdna_wgs.nf path"
  TDNA_NF_PATH="/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/tdna/tdna_wgs.nf"
else
  TDNA_NF_PATH="$DIRname/dna/TDNA/tdna.nf"
fi

cd $workflow_results/FqCramModule/output/cram_output/
mkdir -p $workdir/fq
for j in `ls $workflow_results/FqCramModule/output/cram_output/*fastq_output/*_clean_1.fq.gz |sort -u`;do cp $j $workdir/fq/; done
for j in `ls $workflow_results/FqCramModule/output/cram_output/*fastq_output/*_clean_2.fq.gz |sort -u`;do cp $j $workdir/fq/; done
for j in `ls $workflow_results/FqCramModule/output/cram_output/*fastq_output/*_clean_2.fq.gz |cut -d "/" -f11 |sort -u | sed s/_clean_2.fq.gz//g`;do echo $j; done > $workdir/col3
cd $workdir/fq
for j in `ls *_clean_1.fq.gz |sort -u`;do echo $workdir/fq/$j; done > $workdir/col1
for j in `ls *_clean_2.fq.gz |sort -u`;do echo $workdir/fq/$j; done > $workdir/col2
paste $workdir/col1 $workdir/col2 $workdir/col3 > $workdir/fq.list 
rm -rf $workdir/col1 
rm -rf $workdir/col2 
rm -rf $workdir/col3
cd $workdir

mkdir -p $workdir/output

echo "/mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/bin/nextflow-1 run -bg $TDNA_NF_PATH --fq_list $workdir/fq.list --ref $2 --insert $3 --outdir $workdir/output"
/mnt/lustre/users/sanger-dev/app/bioinfo/dna/env/bin/nextflow-1 run -bg $TDNA_NF_PATH --fq_list $workdir/fq.list --ref $2 --insert $3 --outdir $workdir/output