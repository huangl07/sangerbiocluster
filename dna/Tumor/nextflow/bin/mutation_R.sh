set +e
Rscript --slave $1/bin/panorama_of_genomic_mutations.R --maf all.maf --out ./
