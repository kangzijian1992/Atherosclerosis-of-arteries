rm(list=ls())

library(GEOquery)
library(dplyr)
library(Biobase)
library(tidyr)

####validation the expression in independent dataset

validation<-read.table("Group_Probe_Profile_RAW.txt",header = T,sep = "\t",row.names = 1)
validation[1:10,1:10]

gene_inf<-read.delim("A-MEXP-931.adf.txt",sep = "\t")

head(gene_inf)

gene_inf

rownames(gene_inf)<-gene_inf$Reporter.Name
gene_inf<-gene_inf[rownames(validation),]

table(rownames(validation)==gene_inf$Reporter.Name)

validation<-as.data.frame(validation)
validation$Gene<-gene_inf$Reporter.Database.Entry.hugo.

validation<-aggregate(validation,by = list(validation$Gene),max)
validation[1:10,1:10]

rownames(validation)<-validation$Group.1
validation<-validation[,-1]
validation<-validation[,!colnames(validation)%in%"Gene"]
validation[1:10,1:10]

sample<-read.table("E-MTAB-2055.sdrf.txt",sep = "\t",header = T)

sample

validation[1:10,1:10]
validation<-validation[,grepl(x=colnames(validation),pattern = ".AVG_Signal")]

dim(validation)

colnames(validation)<-sample$Source.Name
sample$Source.Name


sample<-sample[,c(1,24)]
sample

colnames(validation)<-sample$Source.Name
sample<-sample[-27,]

validation<-validation[,sample$Source.Name]

expression<-t(validation)
pdata<-sample
pdata$group<-pdata$Factor.Value.plaque.stability.
validation[1:10,1:10]
library(xlsx)
key_genes<-read.xlsx("keygenes.xlsx",sheetIndex = 1)
key_genes$x

gene_list<-key_genes$x ##set the gene you are interested in
gene_list
i=1
library(ggpubr)  
load("deg.Rdata")
myresults[key_genes$x,]
table(colnames(expression)%in%gene)
pdata$group
gene_list<-gene_list[gene_list%in%colnames(expression)]
gene_list
pdata$group
gene="IGFBP6"
### get the expression data for interested genes
library(stringr)
for (i in 1:length(gene_list)){
  gene<-gene_list[i]
  gene
  
  gene_expression<-expression[,colnames(expression)%in%gene]
  
  gene_expression[1]
  a<-str_count(gene_expression, pattern="[.]")
  gene_expression<-gsub("[.]","",gene_expression)
  gene_expression
  gene_expression<-as.numeric(gene_expression)/(1000^a)
  gene_expression<-log2(gene_expression+1)
  data<-data.frame(expression=gene_expression,group=pdata$group,gene=gene)

  data$group<-factor(data$group,levels = c("stable","ruptured"))

###create box plot
tiff(file=paste0("/Users/yaoyuelin/Desktop/genes_expression/",gene,".tiff"),width = 1200,height = 1500,res = 300)

p<-ggplot(data, aes(x=group, y=expression,fill=group))+geom_boxplot(position=position_dodge(1),outlier.shape = NA)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size =14 ,colour = "black"),panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"),axis.title =element_text(size =25,colour = "black" ),axis.text.y = element_text(size =15,colour = "black" ),legend.text = element_text(size =15,colour = "black" ) ,legend.title =element_text(size =20,colour = "black" ),panel.background = element_rect(fill='transparent', color=NA) )+ylab(expression(paste(log[2],italic("Exp"))))+xlab(" ")
p<-p+facet_wrap( ~gene,scales = "free_y",ncol=1)+theme(
  strip.background = element_rect(color = NA,fill=NA,
                                  size=1),strip.text = element_text(face="bold.italic",size =18),legend.position="none")+stat_compare_means(aes(group =  group),
                     method = "wilcox.test",hide.ns = T)

print(p)
dev.off()


}

myresults["COX7A1",]
myresults["CD68",]



validation[1:10,1:10]

i=1
for (i in 1:length(rownames(validation)))
{
  a<-str_count(validation[i,], pattern="[.]")
  a  
  validation[i,]<-gsub("[.]","",validation[i,])
  validation[i,]<-as.numeric(validation[i,])/(1000^a)
  
  
}
validation[1:10,1:10]




#save(validation,pdata,file="validation_set.Rdata")


