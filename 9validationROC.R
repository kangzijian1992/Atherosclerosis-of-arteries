load("deg.Rdata")
load("remove_batch_effect.Rdata")
load("validation_set.Rdata")
library(glmnet)
dim(combat_mydata)
combat_mydata
dim(Meta_data)

Meta_data$type<-c(rep(c("Stable plaque","Ruptured Plaque"),32),rep("Ruptured Plaque",9),rep("Stable plaque",9),rep("Stable plaque",16),rep("Ruptured Plaque",27),rep("Ruptured Plaque",5),rep("Stable plaque",6),rep(c("Stable plaque","Ruptured Plaque"),4))
table(colnames(combat_mydata)==Meta_data$geo_accession)
deg_exp<-combat_mydata[rownames(myresults)[!myresults$change=="NOT"],]

exp_set<-t(combat_mydata)
combat_mydata[1:10,1:10]
key_genes$x
i="TSPAN6"
library(pROC)
for (i in key_genes$x){

print(i)
gene_exp<-as.data.frame(exp_set[,i])
gene_exp$Status<-Meta_data$type
gene_exp
gene_exp$Status[gene_exp$Status=="Ruptured Plaque"]<-1
gene_exp$Status[gene_exp$Status=="Stable plaque"]<-0
gene_exp$Status<-as.numeric(gene_exp$Status)
gene_exp$Status
colnames(gene_exp)[1]<-"Gene"
gene_exp$Gene<-as.numeric(gene_exp$Gene)
model1 <- glm(Status~Gene, data=gene_exp, family = "binomial")
print(summary(model1))

gene_exp$Base.Probability <- predict(model1, gene_exp, type="response")
gene_exp$Predicted.s <- 1*(gene_exp$Base.Probability>=0.95)
#confusionMatrix(table(gene_exp$Predicted.s,gene_exp$Status),positive = "1")

pdf(file=paste0("/Users/yaoyuelin/Desktop/test_set/",i,".pdf"),width = 3,height = 3)
myroc <- roc(gene_exp$Status, gene_exp$Base.Probability)
plot(myroc,legacy.axes=TRUE,main=i)
text(0.5, 0.8, labels = sprintf(" AUC = %.2f",myroc$auc))
dev.off()

}

key_genes$x
i="SCG3"

load("validation_set.Rdata")

exp_set<-t(combat_mydata)
pdata$group
expression[1:10,1:10]

gene_list<-key_genes$x
gene_list<-gene_list[gene_list%in%rownames(validation)]

for (i in gene_list){
  
  print(i)
  gene_exp<-as.data.frame(exp_set[,i])
  gene_exp$Status<-Meta_data$type
  gene_exp$Status[gene_exp$Status=="Ruptured Plaque"]<-1
  gene_exp$Status[gene_exp$Status=="Stable plaque"]<-0
  gene_exp$Status<-as.numeric(gene_exp$Status)
  gene_exp$Status
  colnames(gene_exp)[1]<-"Gene"
  gene_exp$Gene<-as.numeric(gene_exp$Gene)
  gene_exp
  model1 <- glm(Status~Gene, data=gene_exp, family = "binomial")
  print(summary(model1))
  
  validation
  valiset<-t(validation[i,])
  valiset
  valiset<-as.data.frame(valiset)
  valiset$Status
  valiset$Status<-pdata$group
  
  valiset$Status[valiset$Status=="ruptured"]<-1
  valiset$Status[valiset$Status=="stable"]<-0
  valiset$Status<-as.numeric(valiset$Status)
  colnames(valiset)[1]<-"Gene"
  valiset$Gene<-as.numeric(valiset$Gene)
  
  valiset$Base.Probability <- predict(model1, valiset, type="response")
  valiset$Predicted.s <- 1*(valiset$Base.Probability>=0.8)
  valiset


  pdf(file=paste0("/Users/yaoyuelin/Desktop/test_set/vali_set_",i,".pdf"),width = 3,height = 3)
  myroc <- roc(valiset$Status, valiset$Base.Probability)
  plot(myroc,legacy.axes=TRUE,main=i)
  text(0.5, 0.8, labels = sprintf(" AUC = %.2f",myroc$auc))
  dev.off()
  

  
}




