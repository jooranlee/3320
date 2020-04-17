#you can check your current working directory by typing "getwd()"
setwd("~/Desktop/T5")

#explore the expression values of different genes
geneExpLungbyName  = read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.lung.protein_coding.txt.gz", 
                                row.names = 1,check.names = F)
geneExpLungbyName_t = as.data.frame(t(geneExpLungbyName))

#take logarithm and +1 to get rid of negative no.
log2_gene= log2(geneExpLungbyName_t+1)

library("Hmisc")
gene.rcorr<-rcorr(as.matrix(log2_gene))
#rcorr computes a matrix of pearson's r or spearman's rho rank correlation coefficients
#for all possible pairs of columns of a matrix
gene.p<-gene.rcorr$P
gene.p[upper.tri(gene.p)]<-40

library(reshape2)
gene.pvalues<-melt(gene.p)
head(gene.pvalues)
dim(gene.pvalues)
#filter
library(dplyr)
#eliminate ACE2/TMPRSS2 genes only
#eliminate self-correlating genes
corr_genes<- filter(gene.pvalues, Var2 == "ACE2" | Var2 == "TMPRSS2") %>% filter(Var1 != Var2)
dim(corr_genes)

#filter out the upper matrix values
corr_genes<-filter(corr_genes, value !=40)
dim(corr_genes)
summary(corr_genes)

filtered_genes<- filter(corr_genes, value<0.01)
dim(filtered_genes)
head(filtered_genes)

filtered_ace2<- as.matrix(filter(filtered_genes, Var2 == "ACE2"))
filtered_tmprss2<- as.matrix(filter(filtered_genes, Var2 == "TMPRSS2"))
colnames(filtered_ace2)<- c("Var1", "Var2", "p-value")
colnames(filtered_tmprss2)<- c("Var1", "Var2", "p-value")
head(filtered_ace2)
head(filtered_tmprss2)

ace2<- c(filtered_ace2[,1])
head(ace2)
tmprss2 <- c(filtered_tmprss2[,1])
head(tmprss2)
common<- intersect(ace2, tmprss2)
View(common)

