# RNA-seq
A pipeline records the scripts and command lines for analyzing RNA-seq data for private.

# 1. Trim reads to get high quality clean reads   
`java -Xmx30g -jar trimmomatic-0.36.jar PE -threads 10 sample1_1.fq.gz  sample1_2.fq.gz  sample1_1.clean.fq.gz sample1_1.unpaired.fq.gz sample1_2.clean.fq.gz sample1_2.unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 TOPHRED33 > sample1.log`   
# 2. Build the aligment index using STAR       
`STAR --runMode genomeGenerate --genomeDir /Camel/genome_data --genomeFastaFiles GCF_000767855.1_Ca_bactrianus_MBC_1.0_genomic.fna --runThreadN 24 --sjdbGTFfile GCF_000767855.1_Ca_bactrianus_MBC_1.0_genomic.gtf --sjdbOverhang 149`  
# 3. Alignment    
`#!/bin/sh    
export INDEX=/Camel/genome_data/      
export IN=/Camel/Camel-1/data_release/my_clean      
export OUT=/Camel/bam     
export BAM=//Camel/bam      
export TMP=/Camel/tmp     
for i in \    
1725_HFKKJALXX_L5      
do      
STAR --genomeDir $INDEX --readFilesIn $IN/${i}_1.clean.fq.gz $IN/${i}_2.clean.fq.gz --readFilesCommand zcat --runThreadN 8 --outFilterMultimapNmax 1 --outFilterIntronMotifs RemoveNoncanonical Unannotated --outFilterMismatchNmax 10 --outSAMstrandField intronMotif --outSJfilterReads Unique --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outFileNamePrefix $OUT/${i}     
java -Djava.io.tmpdir=$TMP -jar picard.jar CleanSam I=$OUT/${i}Aligned.out.bam O=$OUT/${i}.STAR.bam     
done`   

