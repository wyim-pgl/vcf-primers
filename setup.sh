
# WS245
mkdir data/WS245
cd data/WS245
curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS245.genomic.fa.gz > WS245/WS245.fa.gz
gunzip WS245.fa.gz
bgzip WS245.fa
samtools faidx WS245.fa.gz