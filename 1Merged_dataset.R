rm(list=ls())

library(GEOquery)
library(dplyr)
library(Biobase)
library(tidyr)
library(xlsx)

"GSE120521"

i="GSE120521"

exprs<-read.table("GSE120521_Athero_RNAseq_FPKM.txt",header = T,row.names = 1)
exprs
gse=getGEO(i,GSEMatrix = TRUE,destdir = ".")

pdata<-pData(gse[[1]])
pdata<-pdata[,c("title","geo_accession")]
pdata
pdata$GEO<-i
exprs
exprs<-exprs[,c(1,7:14)]
exprs<-exprs[!duplicated(exprs$name),]
rownames(exprs)<-exprs$name
exprs<-exprs[,-1]
exprs
pdata$geo_accession
colnames(exprs)<-pdata$geo_accession
pdata$geo_accession

exprs<-as.data.frame(exprs)
exprs

pdata$geo_accession


expression<-exprs
expression[1:5,1:5]
expression$Gene<-rownames(expression)


Expression<-expression
Meta_data<-pdata

#"GSE125771"
#""
i="GSE41571"
for (i in c("GSE41571","GSE163154","GSE111782","GSE43292")) {
  
  print(i)
  gse=getGEO(i,GSEMatrix = TRUE,destdir = ".",getGPL = T,AnnotGPL = T)
  fdata<-fData(gse[[1]])
  head(fdata)
  print(colnames(fdata))
  fdata<-fdata[,c("ID","Gene title","Gene symbol","Gene ID")]
  fdata
  
  #fdata %>% 
  # mutate(`Gene symbol` = strsplit(as.character(`Gene symbol`), "///")) %>% 
  #unnest(`Gene symbol`)
  
  exprs<-exprs(gse[[1]])
  exprs
  pdata<-pData(gse[[1]])
  pdata<-pdata[,c("title","geo_accession")]
  pdata
  pdata$GEO<-i
  exprs<-exprs[rownames(fdata),]
  print(table(rownames(exprs)==rownames(fdata)))
  exprs<-as.data.frame(exprs)
  exprs$Gene<-fdata$`Gene symbol`
  
  
  
  exprs<-exprs %>% 
    mutate(`Gene` = strsplit(as.character(`Gene`), "///")) %>% 
    unnest(`Gene`)
  
  exprs<-data.frame(exprs)
  expression<-aggregate(exprs,by = list(exprs$Gene),max)
  

  rownames(expression)<-expression$Group.1
  expression<-expression[,-1]
  
  print(dim(expression))
  
  Expression<-inner_join(Expression,expression,by="Gene")
  
  Meta_data<-rbind(pdata,Meta_data)
}

Meta_data$title

Meta_data$geo_accession

Expression[1:10,1:10]

Expression[1:10,"Gene"]
rownames(Expression)<-Expression$Gene
Expression<-Expression[,!colnames(Expression)%in%"Gene"]
Expression<-Expression[,Meta_data$geo_accession]
save(Meta_data,Expression,file="Merge_expression.Rdata")

