rm(list=ls())

library(GEOquery)
library(dplyr)
library(Biobase)
library(tidyr)

load("Merge_expression.Rdata")

Expression[1:10,1:10]

Expression<-as.matrix(Expression)

Meta_data$GEO
Meta_data$title

Expression<-Expression[,Meta_data$geo_accession]
unique(Meta_data$GEO)
Expression[1:10,1:10]
number<-table(Meta_data$GEO)[unique(Meta_data$GEO)]
number
colours <- c(rep("yellow",number[1]),rep("red",number[2]),rep("blue",number[3]),rep("purple",number[4]),rep("orange",number[5]))
names(colours)<-Meta_data$GEO

###boxplot before 
Meta_data$GEO
dim(Expression)
pdf("boxplot_before.pdf")
boxplot(Expression,col=colours,xaxt="n",outline=F,frame = FALSE,cex.axis=1.5)
legend("topright",inset = 0.01, legend=unique(Meta_data$GEO),
       col=unique(colours),lty=1, cex=1,pch = 15,box.col = "white",bg = "white")#bty = "n"
title(xlab="Sample",line=2.5,ylab="Gene Expression",cex.lab = 2)
dev.off()

###pca before
library(factoextra)
library(FactoMineR)
dat.pca<-PCA(t(Expression), graph = FALSE)
Meta_data$GEO

pdf("pca_before.pdf",width = 5,height = 4)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Meta_data$GEO, # color by groups
                         palette = unique(colours[order(Meta_data$GEO)]),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups")
print(pca_plot)
dev.off()

#######



Meta_data$title
Meta_data$GEO
Meta_data$type<-Meta_data$title

####remove batch effect
library("sva")
#BiocManager::install("sva")
Expression[1:10,1:10]

exp<-Expression
pheno<-Meta_data[,c("type","GEO")]

table(Meta_data$GEO)
batch<-Meta_data$GEO
modcombat<-model.matrix(~1, data=pheno)
combat_mydata= ComBat(dat=exp, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


#save(combat_mydata,Meta_data,file="remove_batch_effect.Rdata")

###boxplot after
pdf("boxplot_after.pdf")
boxplot(combat_mydata,col=colours,xaxt="n",outline=F,frame = FALSE,cex.axis=1.5)
legend("topright",inset = 0.01, legend=unique(Meta_data$GEO),
       col=unique(colours),lty=1, cex=1,pch = 15,box.col = "white",bg = "white")#bty = "n"
title(xlab="Sample",line=2.5,ylab="Gene Expression",cex.lab = 2)
dev.off()
###pca
library(factoextra)
library(FactoMineR)
dat.pca<-PCA(t(combat_mydata), graph = FALSE)

pdf("pca_after.pdf",width = 5,height = 4)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = factor(Meta_data$GEO), # color by groups
                         palette = unique(colours[order(Meta_data$GEO)]),
                         addEllipses = F, # Concentration ellipses
                         legend.title = "Groups")


print(pca_plot)
dev.off()


