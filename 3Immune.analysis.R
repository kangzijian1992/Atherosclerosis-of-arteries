load("remove_batch_effect.Rdata")

######write table for input
combat_mydata<-data.frame(combat_mydata)
ematrix<-combat_mydata
cols =  c("GeneSymbol",colnames(ematrix))
rows <- rownames(ematrix)
ematrix <- cbind(rows,ematrix)
ematrix <- ematrix[which(ematrix[,1] != "NA"),] #remove NAs
ematrix <- ematrix[order(ematrix[,1]),] #sort by gene name 
ematrix <- rbind(cols, ematrix)
ematrix[1:10,1:10]
write.table(ematrix,file="NormalizedExpressionArray.customCDF.txt",sep="\t", col.names=F, row.names=F,quote=FALSE)


##### get result from website: CIBERSORT https://cibersort.stanford.edu/
immu.result<-as.data.frame(read.csv("CIBERSORT.Output_Job5.csv",header = T,row.names = 1))[,1:22]
immu.result

sum<-colSums(immu.result)
sum

immu.result<-immu.result[,!sum==0,]
dim(immu.result)
#immu.result$Group<-Meta_data$type
immu.result$Sample<-Meta_data$geo_accession
table(rownames(immu.result)==Meta_data$geo_accession)

####

library(dplyr)
library(reshape2)
library(ggpubr)
TME_New = melt(immu.result)
head(TME_New)

immu.result
Meta_data$type<-c(rep(c("Stable plaque","Ruptured Plaque"),32),rep("Ruptured Plaque",9),rep("Stable plaque",9),rep("Stable plaque",16),rep("Ruptured Plaque",27),rep("Ruptured Plaque",5),rep("Stable plaque",6),rep(c("Stable plaque","Ruptured Plaque"),4))
table(Meta_data$type)
TME_New$Group=Meta_data$type
TME_New
colnames(TME_New)=c("Sample","Celltype","Composition","Group")
TME_New$Group

plot_order = TME_New[TME_New$Group=="Ruptured Plaque",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)


if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12))}


box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",hide.ns = T)

box_TME
box_TME;ggsave("TME.pdf",box_TME,height=12,width=20,unit="cm")

TME_New
result<-compare_means(Composition~Group,data =TME_New,group.by="Celltype",method = "wilcox.test" )
result<-as.data.frame(result)

significant<-result[result$p.signif!="ns",] 
significant
box_TME

head(TME_New)
TME_New<-TME_New[order(TME_New$Group),]



library(tidyverse)

p<-ggplot(TME_New, mapping = aes(x =Sample , y = Composition, fill =  Celltype)) + scale_y_continuous(labels = scales::percent)+
  geom_bar(position= "stack",stat = "identity")+theme_bw()+theme(axis.text.x = element_blank(),panel.spacing = unit(0.1, "lines"))
p<-p+facet_grid(. ~ Group,space = "free",scales = "free_x")+theme(strip.text.x = element_text(
  size = 10,face = "bold.italic"),strip.background = element_rect(colour=NA, fill="white",size=15, linetype="solid"))+ylab("Relative Percent")
p
ggsave("composition.pdf",p,height=10,width=22,unit="cm")
#####
library(pheatmap)
head(immu.result)


annotation=data.frame(Group=Meta_data$type)
rownames(annotation)=rownames(immu.result)
table(rownames(immu.result)==Meta_data$geo_accession)



ann_colors = list(
  Group = c("red", "blue")
)


pdf("heatmap.pdf",width = 8,height = 5)
pheatmap(t(immu.result[Meta_data$geo_accession[order(Meta_data$type)],1:22]),annotation_col = annotation,cluster_cols = F,show_colnames = F)
dev.off()
Meta_data$type
immu.result
mat<-immu.result[Meta_data$geo_accession[Meta_data$type=="Ruptured Plaque"],1:22]
mat

mat_cor<-cor(mat,method = "spearman")
mat_cor

library(corrplot)
pdf("corplot.pdf",width = 15,height = 12)
corrplot(mat_cor,method = 'number',order = 'hclust',tl.cex = 1.5)
dev.off()
######## step 1 Build the design matrix
library(factoextra)
library(FactoMineR)
dat.pca<-PCA(t(combat_mydata), graph = FALSE)
dat.pca
pdf("pca_group.pdf",width = 4,height = 3)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = factor(Meta_data$type), # color by groups
                         palette = c("blue","red"),
                         addEllipses = F, # Concentration ellipses
                         legend.title = "Groups")


print(pca_plot)
dev.off()


pdf("pca_immu.pdf",width = 4,height = 3)
dat.pca<-PCA(immu.result[1:21], graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = factor(Meta_data$type), # color by groups
                         palette = c("blue","red"),
                         addEllipses = F, # Concentration ellipses
                         legend.title = "Groups")


print(pca_plot)
dev.off()


######laosso
#install.packages("glmnet")
library(glmnet)
immu.result[1:10,1:10]
x=data.frame(immu.result[,1:22])
y=Meta_data$type
(immu.result)
table(rownames(immu.result)==Meta_data$geo_accession)
Meta_data$type
set.seed(1000)
cvfit<-cv.glmnet(as.matrix(x),y,family = "binomial")
cvfit$lambda.min
cvfit$lambda.1se
pdf("cvfit.pdf",width = 5,height = 4)
plot(cvfit)
text(-4,1.3,"lambda.min: 0.02 \nlambda.1se: 0.05")
dev.off()


#model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cvfit$lambda.min)
fit<-glmnet(x,y,family = "binomial",alpha = 1,lambda = cvfit$lambda.min)
lasso_select<-rownames(fit$beta)[as.numeric(fit$beta)!=0]
significant$Celltype
lasso_select

fit3<-cvfit

fit3.coef.lambda.min<-coef(fit3,s=fit3$lambda.min)
fit3.min.out<-fit3.coef.lambda.min[which(fit3.coef.lambda.min!=0),]
fit3.min.out
fit3.min.out<-round(fit3.min.out,4)
write.xlsx(fit3.min.out,file="lasso_immune_cell.coef.min.xlsx")




key_immune<-intersect(lasso_select,significant$Celltype)
key_immune
fit<-glmnet(x,y,alpha = 1,nlambda = 1000,family = "binomial")
pdf("fit.pdf",width = 5,height = 4)
plot(fit,xvar="lambda",label=F)
dev.off()

y

library(randomForest)
y[y=="Ruptured Plaque"]<-1
y[y=="Stable plaque"]<-0
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
rf_output=randomForest(x=x, y=factor(y),importance = TRUE,mtry = best_mtry, ntree = 10000, proximity=TRUE )
y
model<-rf_output
model$err.rate
oob.error.data <- data.frame(
  Trees=rep(1:nrow(model$err.rate), times=3),
  Type=rep(c("OOB", "Ruptured Plaque", "Stable plaque"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"],
          model$err.rate[,"Stable plaque"],
          model$err.rate[,"Ruptured Plaque"]))
pdf("tree.pdf",width = 5,height = 3)
p<-ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +theme_bw()+
  geom_line(aes(color=Type))
print(p)
dev.off()


#ntree=400

rf_output=randomForest(x=x, y=factor(y),importance = TRUE,mtry = best_mtry, ntree = 7500, proximity=TRUE )
rf_output$err.rate
importance(x=rf_output)
rf_output

rf_importances=importance(x=rf_output)
rf_output
choose_gene1=rownames(tail(rf_importances[order(rf_importances[,3]),],10))
choose_gene2=rownames(tail(rf_importances[order(rf_importances[,4]),],10))

key_immu<-intersect(unique(c(choose_gene1,choose_gene2,lasso_select)),significant$Celltype)
choose_gene2
key_immu
significant$Celltype
#save(key_immu,file="key_immu.Rdata")

pdf("important.pdf",width = 10)
varImpPlot(rf_output)
dev.off()





