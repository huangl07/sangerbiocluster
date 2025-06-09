 plink --vcf 1000.vcf  --recode structure  --allow-extra-chr --biallelic-only
 less -S plink.recode.strct_in|perl -ne '$n++;print $_ if($n > 2);' > pop.structure
 N=`wc -l pop.structure`
 L=`less -S pop.structure|head -n 1|awk '{print NF-2}'`


for i in {1..10}
do
    for j in {1..10}
    do
        echo -e  "structure -m  ~/sg-users/long.huang/STRUCTURE/mainparams -e ~/sg-users/long.huang/STRUCTURE/extraparams -K ${i} -i pop.structure  -N ${N} -L ${L} -o K${i}_run${j}" 
    done
done|xargs -i -P 8 bash -c '{}'


structureHarvester.py --dir=./ --out=./ --clumpp

for i in {1..10}
do
    echo -e "CLUMPP -i K${i}.indfile -o K${i}.outfile -j K${i}.miscfile -p K${i}.popfile -R 10 -K ${i} -c ${N}"
done|xargs -i -P 8 bash -c '{}'

cp  ~/sg-users/long.huang/STRUCTURE/clumpp.para paramfile
for i in {1..10}
do
    echo -e "less -S K${i}.outfile |perl -ne '\$n++;@a=split;print join("\t",\$n.\$a[4],@a[5..\$#a]),"\t1\n"' > K${i}.newpop"
done|xargs -i -P 8 bash -c '{}'

