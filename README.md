# RNA-seq
A pipeline records the scripts and command lines for analyzing RNA-seq data for private.

# 1. Trim reads to get high quality clean reads   
`java -Xmx30g -jar trimmomatic-0.36.jar PE -threads 10 sample1_1.fq.gz  sample1_2.fq.gz  sample1_1.clean.fq.gz sample1_1.unpaired.fq.gz sample1_2.clean.fq.gz sample1_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 > sample1.log`   
# 2. Build the aligment index using STAR    
`STAR --runMode genomeGenerate --genomeDir /Camel/genome_data --genomeFastaFiles /Camel/genome_data/GCF_000767855.1_Ca_bactrianus_MBC_1.0_genomic.fna --runThreadN 24 --sjdbGTFfile /Camel/genome_data/GCF_000767855.1_Ca_bactrianus_MBC_1.0_genomic.gtf --sjdbOverhang 149`    

