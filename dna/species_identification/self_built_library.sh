# $1 需要创建的数据库名字
# $2 需要创建的数据库包含的物种拉丁文list，一行一个


source ~/sg-users/yiwei.tang/tyw.rc
taxonkit name2taxid  -o name2taxid -j 8 $2
cut -f 2 name2taxid |sed '/^$/d' >taxid.list
cat /mnt/lustre/users/sanger-dev/sg-users/yiwei.tang/test_data/nt/nucl.accession2taxid | csvtk -t grep -f taxid -P taxid.list | csvtk -t cut -f accession.version > taxa.giid.acc.txt
blastdbcmd -db /mnt/lustre/users/sanger-dev/app/database/nt/nt/nt_20230606/nt -entry_batch taxa.giid.acc.txt -out $1.fa

export PATH=$PATH:/mnt/lustre/users/sanger-dev/app/bioinfo/dna/kraken2-master/kraken2-master/
/mnt/lustre/users/sanger-dev/app/bioinfo/dna/kraken2-master/kraken2-master/scripts/kraken2-build --add-to-library $1.fa --threads 8 --db $1
cp -rl /mnt/lustre/users/sanger-dev/app/database/nt/nt/nt_20230608_Kraken/taxonomy ./$1
/mnt/lustre/users/sanger-dev/app/bioinfo/dna/kraken2-master/kraken2-master/scripts/kraken2-build --build --threads 50 --db $1
/mnt/lustre/users/sanger-dev/app/bioinfo/dna/Bracken/bracken-build -d $1 -t 8 -k 35 -l 150

