# PASA基于组装好的unigene注释基因组

```bash
# 准备gene集，核酸序列，unigene.fa和参考基因组序列ref.fa，运行PASA
nextflow run -bg pasa-dsl2.nf --ref ref.fa --transcript Trinity.fasta --outdir result
```

# Braker3注释基因组

```bash
# 无转录组数据
nextflow run -bg braker3-dsl2.nf --ref ref.fa --db /mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Arthropoda.fa --outdir result


# 提供转录组数据
nextflow run -bg braker3-dsl2.nf --ref ref.fa --db /mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Arthropoda.fa --fqlist fq.list --outdir result
```

> **注1**
>
> --fqlist
> fq.list 为三列，第一列样品名，第二列reads1，第二列reads2。
>
> **注2**
>
> --db    <file>  database, 目前提供5个，按需选择
>
> 1. 节肢动物门：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Arthropoda.fa
>
> 2. 脊椎动物门：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Vertebrata.fa
>
> 3. 植物界：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Viridiplantae.fa
>
> 4. 真菌：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Fungi.fa
>
> 5. 真核生物：/mnt/lustre/users/sanger-dev/sg-users/pengyu.guo/docker_images/orthodb/Eukaryota.fa
