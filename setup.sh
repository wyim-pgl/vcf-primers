
# WS245
mkdir reference
mkdir reference/WS245
cd reference/WS245
curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/genomic/c_elegans.PRJNA13758.WS245.genomic.fa.gz > WS245.fa.gz
gunzip -f WS245.fa.gz
bgzip --stdout WS245.fa > WS245.fa.gz

# Make bwa index
bwa index WS245.fa.gz

# Make samtools db.
samtools faidx WS245.fa.gz

# Generate blast db.
makeblastdb -in WS245.fa -dbtype=nucl 