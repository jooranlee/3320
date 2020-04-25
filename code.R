#you can check your current working directory by typing "getwd()"
setwd("~/Desktop/T5")

donorInfo = read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.lung.protein_coding.txt", sep = "\t",
                       header = T, stringsAsFactors = F)
#explore the expression values of different genes
geneExpLungbyName  = read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.lung.protein_coding.txt.gz", 
                                row.names = 1,check.names = F)
geneExpLungbyName_t = as.data.frame(t(geneExpLungbyName))

#finding correlated genes to ACE2
mylist1<-vector()
numberofgenes1<-0

for (i in 1:ncol(geneExpLungbyName_t)){
  x=donorInfo$GeneName[i]
  if(x=="ACE2"){
    next
  }
  ctp1S=cor.test(x=geneExpLungbyName_t$ACE2, y=geneExpLungbyName_t[[x]],method = "spearman" )
  ctp1P=cor.test(x=geneExpLungbyName_t$ACE2, y=geneExpLungbyName_t[[x]],method = "pearson" )
  pvalue1P=ctp1P$p.value
  pvalue1S=ctp1S$p.value
  if (is.na(pvalue1P)||is.na(pvalue1S)){
    next
  }
  if(pvalue1P<5.189144e-07 && pvalue1S<5.189144e-07){
    mylist1<-append(mylist1,x)
    numberofgenes1=numberofgenes1+1
  }
} 
#check View(mylist1)

#find the correlated genes to TMPRSS2
mylist2<-vector()
numberofgenes2<-0

for (i in 1:ncol(geneExpLungbyName_t)){
  x=donorInfo$GeneName[i]
  if(x=="TMPRSS2"){
    next
  }
  ctpl2S=cor.test(x=geneExpLungbyName_t$TMPRSS2, y=geneExpLungbyName_t[[x]],method = "spearman" )
  ctpl2P=cor.test(x=geneExpLungbyName_t$TMPRSS2, y=geneExpLungbyName_t[[x]],method = "pearson" )
  pvalue2P=ctpl2P$p.value
  pvalue2S=ctpl2S$p.value
  if (is.na(pvalue2P)||is.na(pvalue2S)){
    next
  }
  if(pvalue2P<0.001/19271 && pvalue2S<0.001/19271){
    mylist2=append(mylist2,x)
    numberofgenes2=numberofgenes2+1
  }
}   

# check View(mylist2)
#get common gene list that correlates to both ACE2 and TMPRSS2 (intersect)
cat(sapply(mylist3, toString), file="genelist.txt", append=FALSE, sep="\n")
mylist3=intersect(mylist1, mylist2)
