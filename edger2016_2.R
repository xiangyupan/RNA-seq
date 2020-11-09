rm(list=ls())
library(edgeR)
#Data import:
setwd("C:/Users/Jack Zhou/Desktop/duck_code/MITF")
xpre <- read.delim("wild_Z2_F1_NCBI2016.readscount.txt",row.names="Gene")
xpre[1:10,]
x <- xpre[,c(22,23,24,25,26,27)]
#x <- xpre
#[,c(7,8,9,10,11)]
dim(x)
head(x)
#group <- factor(c(1,1,1,2,2,2,1,1,2,2,2,1,2,1,2,1,2,1,2,1,2,1,1,2,2,2,1,1,2,2,2,1,2,1,2,1,2,1,2,1,2,1,1,2,2,2,1,1,2,2,2,1,2,1,2,1,2,1,2))
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
y$sample
head(y$counts)
summary(y$counts)
#Reset the library sizes after filter;
y$samples$lib.size <- colSums(y$counts)
y$samples
#Normalization;
y <- calcNormFactors(y)
    y$pseudo.counts[1:10,]
#write.table(y$pseudo.counts, 'Z2_wild_F1.norm2016.readscount.txt', sep = '\t', row.names = T, col.names = T, quote = F)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y, trend="none")
plotBCV(y, cex=0.4)
et <- exactTest(y)
topTags(et, n=20)
detags <- rownames(topTags(et, n=20))
cpm(y)[detags,]
#-1, 0 and 1 are for down-regulated, non-differentially expressed and up-regulated
summary(de <- decideTestsDGE(et, p=0.05, adjust="BH"))
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")
#write.table(detags, 'genes_DE_edgeR_bremuscle.txt', sep = '\t', row.names = F, col.names = T, quote = F)
write.table(topTags(et, n=388), 'genes_DE_edgeR_fat.txt', sep = '\t', row.names = T, col.names = T, quote = F)


从这开始跑！！！
##########################
library("limma")
library("edgeR")
a<-read.table("Colon.sika_deer.atleast1.ID.txt",header=T) 
head(a)
y<-DGEList(counts=a[,2:16],genes=a[,1])
#left<-rowSums(cpm(y)>1)>=4
#y<-y[left,]
y<-DGEList(counts=y$counts,genes=y$genes)
y<-calcNormFactors(y)
group<-factor(c(1,1,1,1,1,2,2,2,2,2))
design <- model.matrix(~group)
y <-DGEList(counts=a[,2:16],genes=a[,1])
y<-estimateGLMCommonDisp(y,design,verbose=TRUE)
y<-estimateGLMTrendedDisp(y, design)
y<-estimateGLMTagwiseDisp(y, design)
fit <- glmQLFit(y,design)
####qif是一个很复杂的文件格式
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
fit<-glmFit(y, design)
lrt<-glmLRT(fit)
topTags(lrt)
##火山图
summary(de<-decideTestsDGE(qlf))
detags<-rownames(y)[as.logical(de)]
plotSmear(qlf, de.tags=detags)
abline(h=c(-2,2),col='blue')
###输出table
table<-cbind(qlf$genes,qlf$table)
write.table(table, 'Colon.DEG.txt', sep = '\t', row.names = T, col.names = T, quote = F)

###file-change dir-选择文件夹
a<-read.table("Colon.DEG.txt",header=T)
b <-a$PValue
p.adjust(b,method="fdr",n=26753)
c<-p.adjust(b,method="fdr",n=26753)
table<-(cbind(a,c))
write.table(table,"Colon.DEG.after.fdr.tsv",sep = '\t')