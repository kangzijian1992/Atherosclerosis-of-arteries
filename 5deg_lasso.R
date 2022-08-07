load("deg.Rdata")
table(myresults$change)
load("remove_batch_effect.Rdata")
library(glmnet)
dim(combat_mydata)
combat_mydata
dim(Meta_data)


table(colnames(combat_mydata)==Meta_data$geo_accession)

deg_exp<-combat_mydata[rownames(myresults)[!myresults$change=="NOT"],]


deg_exp[1:10,1:10]
x

Meta_data$type<-c(rep(c("Stable plaque","Ruptured Plaque"),32),rep("Ruptured Plaque",9),rep("Stable plaque",9),rep("Stable plaque",16),rep("Ruptured Plaque",27),rep("Ruptured Plaque",5),rep("Stable plaque",6),rep(c("Stable plaque","Ruptured Plaque"),4))

x=t(deg_exp)
y=Meta_data$type
y
y[Meta_data$type=="Ruptured Plaque"]<-1
y[Meta_data$type=="Stable plaque"]<-0
y<-as.numeric(y)
y
set.seed(1234)
fit3 <- cv.glmnet(x=x, y=factor(y),type.measure="deviance", alpha = 1,family = "binomial")


pdf("cvfit_genes.pdf",width = 5,height = 4)
plot(fit3)
text(-4,1.4,"lambda.min: 0.05 \nlambda.1se: 0.10")
dev.off()

fit3$lambda.min
fit3$lambda.1se

fit3.coef.lambda.1se<-coef(fit3,s=fit3$lambda.1se)
fit3.1se.out<-fit3.coef.lambda.1se[which(fit3.coef.lambda.1se!=0),]
fit3.1se.out
fit3.1se.out<-round(fit3.1se.out,4)

fit3.coef.lambda.min<-coef(fit3,s=fit3$lambda.min)
fit3.min.out<-fit3.coef.lambda.min[which(fit3.coef.lambda.min!=0),]
fit3.min.out
fit3.min.out<-round(fit3.min.out,4)
library(xlsx)
write.xlsx(fit3.min.out,file="lasso.coef.min.xlsx")

fit<-glmnet(x,y,family = "binomial")
pdf("fit_genes.pdf",width = 5,height = 4)
plot(fit,xvar="lambda",label=F)
dev.off()


cfit <- cv.glmnet(x, y, family = "binomial", type.measure = "auc", 
                  keep = TRUE)
rocs <- roc.glmnet(cfit$fit.preval, newy = y)
best <- cvfit$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")
#https://glmnet.stanford.edu/articles/glmnet.html#assessing-models-on-test-data-1

#####
lasso.prob <- predict(fit3, newx=x , s=c(fit3$lambda.min,fit3$lambda.1se),family = "binomial",type="response")
lasso.prob
re=cbind(y ,lasso.prob)
#合并样本预测值与真实值
(re)
re
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)
library(ggpubr) 

p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p1


library(ROCR)
re
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
#求得AUC值
pdf("auo_deg_lasso.pdf",width = 5,height = 5)
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="red") 
#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
dev.off()
# 加AUC值






#####
library(randomForest)

y
y=Meta_data$type
errset=1:(length(y)-1)
#library(e1071) svm
for (i in 1:length(y)-1){
  set.seed(350)
  mtry_fit=randomForest(x=x, y=factor(y),mtry=i)
  err=mean(mtry_fit$err.rate)
  print(err)
  errset[i]=err
}

errset
best_mtry=which(errset==min(errset),arr.ind=TRUE)
best_mtry

rf_output=randomForest(x=x, y=factor(y),importance = TRUE,mtry = best_mtry, ntree = 1500, proximity=TRUE )

model<-rf_output

oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=3),
  Type=rep(c("OOB", "Ruptured Plaque", "Stable plaque"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
          model$err.rate[,"Stable plaque"],
          model$err.rate[,"Ruptured Plaque"]))
pdf("tree_genes.pdf",width = 5,height = 3)
p<-ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +theme_bw()+
  geom_line(aes(color=Type))
print(p)
dev.off()


#ntree=750

rf_output=randomForest(x=x, y=factor(y),importance = TRUE,mtry = best_mtry, ntree = 750, proximity=TRUE )
rf_output$err.rate
importance(x=rf_output)
rf_output

rf_importances=importance(x=rf_output)
rf_output
choose_gene1=rownames(tail(rf_importances[order(rf_importances[,3]),],30))
choose_gene2=rownames(tail(rf_importances[order(rf_importances[,4]),],30))


pdf("important_genes.pdf",width = 10)
varImpPlot(rf_output)
dev.off()



#library(e1071) svm


lass<-names(fit3.min.out)[-1]
choose_gene2
choose_gene1
lass

key_genes<-unique(c(choose_gene2,choose_gene1,lass))
key_genes

write.xlsx(key_genes,file="keygenes.xlsx")

x <- list(
  RF1 = choose_gene2, 
  RF2 = choose_gene1, 
  LASSO= lass
)

library(ggvenn)
pdf("ven.pdf",width = 4,height = 4)
p<-ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
print(p)
dev.off()


#key_genes<-key_genes$x
expre_heatmap<-combat_mydata[key_genes,]

expre_heatmap[1:10,1:10]

library(pheatmap)

table(colnames(expre_heatmap)==Meta_data$geo_accession)

annotation=data.frame(Group=Meta_data$type)
rownames(annotation)=colnames(expre_heatmap)
table(colnames(expre_heatmap)==Meta_data$geo_accession)

annotation

immu.result<-as.data.frame(read.csv("CIBERSORT.Output_Job5.csv",header = T,row.names = 1))[,1:22]
ann_colors = list(
  Group = c("red", "blue")
)
library(RColorBrewer)
breaksList = seq(-1, 1, by = 0.1)
#expre_heatmap<-log2(expre_heatmap+1)
pdf("heatmap_key_genes.pdf",width = 8,height = 8)
pheatmap(expre_heatmap[,order(Meta_data$type)],annotation_col = annotation,cluster_cols = F,cluster_rows = T,show_colnames = F,border_color = NA,scale = 'row',breaks = breaksList,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
dev.off()

load("key_immu.Rdata")
key_immu
Meta_data$type
expre_heatmap<-combat_mydata[key_genes,Meta_data$geo_accession[Meta_data$type=="Ruptured Plaque"]]

#immu_mtr<-immu.result[Meta_data$geo_accession[Meta_data$type=="Meningioma"],key_immu[-7]]
immu_mtr<-immu.result[Meta_data$geo_accession[Meta_data$type=="Ruptured Plaque"],key_immu]

expre_heatmap
immu_mtr

mat<-cor(cbind(immu_mtr,t(expre_heatmap)),method = "spearman")
mat
#save(mat,file="mat_gene_immune.Rdata")
pdf("corplot_gene_cells.pdf",width = 10,height = 10)
corrplot(mat,method = 'number',tl.cex = 0.7,cl.cex = 1.5,number.cex = 0.45)
dev.off()

pdf("corplot_gene_cells1.pdf",width = 10,height = 10)
corrplot(mat,tl.cex = 0.7,cl.cex = 1.5,number.cex = 0.45)
dev.off()



write.xlsx(mat,file="cor_mat_gene_cell.xlsx")

melt<-melt(mat)
melt<-melt[!melt$value==1,]
melt<-melt[order(melt$value,decreasing = T),]

melt<-melt[melt$Var1%in%key_immu,]
tail(melt,20)
head(melt,20)
selected_genes<-c("PCOLCE2","CD68","CNN1")
selected_cells<-c("Plasma.cells","Macrophages.M0","T.cells.CD8")

for (i in 1:3){

genes<-data.frame(cells=immu_mtr[,selected_cells[i]],genes=expre_heatmap[selected_genes[i],])
pdf(file =paste0( "correlation_",selected_genes[i],".pdf"),width = 3.5, height = 3.5)
par(mar=c(5.1,6,6,10))
DTI1 <-ggscatter(genes, x = "cells" , y = "genes", font.label = 14, repel = T,
                 
                 conf.int = T,  , size = 2, cor.coef.size = 4, add = "reg.line",
                 
                 cor.coef = T,cor.method = "spearman", #label = rownames(comb),
                 
                 xlab = paste0(selected_cells[i]) , ylab =(selected_genes[i]), cor.coeff.args = list(label.sep = "\n",col="red"),
                 
                  cex.lab=1,cex.axis = 0.8,)

p<-DTI1+font("xlab", size =20 , color = "black")+
  font("ylab", size = 20, color = "black")+
  font("xy.text", size = 10, color = "black")+theme(legend.position='none')

print(p)
dev.off()
}



mat<-cor(immu_mtr,t(expre_heatmap),method = "spearman")
mat

graph<-melt(mat)
graph$abs_value<-abs(graph$value)
graph<-graph[graph$abs_value>0.1,]
graph
save(graph,file="graph.Rdata")
length(unique(graph$Var2))
degree<-table(graph$Var2)
df<-as.data.frame(degree[order(degree,decreasing = T)])
df
pdf(file =paste0( "degree.pdf"),width = 6, height = 3)
p<-ggplot(data=df, aes(x=Var1, y=Freq,group=1)) +theme_bw()+
  geom_line()+theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 7.0))+ylab("Degree")+xlab("Gene")+
  geom_point()
print(p)
dev.off()
library(igraph)

load("merged.Rdata")
merged<-merged[merged$geneName%in%df$Var1[df$Freq>1],]
dim(merged)
merged$sum<-merged$PITA+merged$RNA22+merged$miRmap+merged$microT+merged$miRanda+merged$PicTar+merged$TargetScan
merged_new<-merged[merged$sum>4,]
dim(merged_new)
write.xlsx(merged_new[,c(2,4)],file="miRNA_gene_interaction.xlsx")
write.xlsx(merged_new,file="miRNA_gene_interaction_all.xlsx")
merged_new[,c(2,4)]
