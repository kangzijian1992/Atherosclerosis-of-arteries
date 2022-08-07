library(limma)
load('remove_batch_effect.Rdata')
combat_mydata[1:10,1:10]

Meta_data$type<-c(rep(c("Stable plaque","Ruptured Plaque"),32),rep("Ruptured Plaque",9),rep("Stable plaque",9),rep("Stable plaque",16),rep("Ruptured Plaque",27),rep("Ruptured Plaque",5),rep("Stable plaque",6),rep(c("Stable plaque","Ruptured Plaque"),4))
Meta_data$title





group_list<-Meta_data$type
Meta_data$type
design <- model.matrix(~0+group_list)
design
colnames(design) <- c("Ruptured_Plaque","Stable_plaque")
rownames(design)=colnames(combat_mydata)

dim(design)


##########step 2Contrasts Matrix
contrastmatrix <- makeContrasts(Ruptured_Plaque-Stable_plaque, levels=design)#This instructs Limma which comparisons to make
contrastmatrix

###############step 3Fit the Linear Model
fit <- lmFit(combat_mydata, design) ##issue these commands to fit the model and make the contrasts
fit2 <- contrasts.fit(fit, contrastmatrix)
fit2 <- eBayes(fit2)
fit2
myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(combat_mydata))
myresults 
logFC_cutoff=1.5
#myresults[order(myresults$logFC,decreasing = T),]
myresults$change = as.factor(ifelse(myresults$adj.P.Val < 0.05 & abs(myresults$logFC) > logFC_cutoff,
                                    ifelse(myresults$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
table(myresults$change)
up_regulated<-myresults[myresults$change=="UP",]
up_regulated<-up_regulated[order(up_regulated$logFC,decreasing = T),]
down_regulated<-myresults[myresults$change=="DOWN",]
down_regulated<-down_regulated[order(down_regulated$logFC,decreasing = F),]
myresults

save(myresults,file="deg.Rdata")

library(xlsx)
write.xlsx(myresults,file="limma_DEG.xlsx")
summary(myresults)
myresults$display<-NA

myresults$display[rownames(myresults)%in%c(rownames(up_regulated)[1:15],rownames(down_regulated)[1:15])]<-rownames(myresults)[rownames(myresults)%in%c(rownames(up_regulated)[1:15],rownames(down_regulated)[1:15])]

library(ggrepel)

pdf(paste0("vocanol.pdf"),width = 10,height = 6)
p<-ggplot(data= myresults, aes(x=logFC, y=-log10(adj.P.Val), col=change,label=display)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +xlim(c(-50,50))+
  scale_color_manual(values=c("blue", "black", "red")) +ylab("-log10(FDR)")+
  geom_vline(xintercept=c(-1.5,1.5), col="red") +theme(legend.title = element_blank(),legend.text = element_text(size=15),axis.title = element_text(size=15))+
  geom_hline(yintercept=-log10(0.05), col="red")
print(p)

dev.off()

c(rownames(up_regulated)[1:10],rownames(down_regulated)[1:10])

expre_heatmap<-combat_mydata[c(rownames(up_regulated)[1:20],rownames(down_regulated)[1:20]),]

expre_heatmap[1:10,1:10]

library(pheatmap)

table(colnames(expre_heatmap)==Meta_data$geo_accession)

annotation=data.frame(Group=Meta_data$type)
rownames(annotation)=colnames(expre_heatmap)
table(colnames(expre_heatmap)==Meta_data$geo_accession)

annotation


ann_colors = list(
  Group = c("red", "blue")
)
breaksList = seq(-1, 1, by = 0.1)
#expre_heatmap<-log2(expre_heatmap+1)
library(RColorBrewer)
pdf("heatmap1.pdf",width = 8,height = 5)
pheatmap(expre_heatmap[,order(Meta_data$type)],annotation_col = annotation,cluster_cols = F,
         cluster_rows = T,show_colnames = F,border_color = NA,scale = 'row',breaks = breaksList,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
dev.off()

#####go and kegg
rm(list=ls())
########################################
###########load packages################
########################################

library(org.Hs.eg.db)
library(GOplot)
library(clusterProfiler)
library("KEGGREST")
library(igraph)
library(ggraph)
library(scales)
library(data.table)
library(network)
library(tidygraph)
library(dplyr)
#install.packages('GOplot')
library(xlsx)
library(biomaRt)
#####################################
########load gene list###############
#####################################
DEG<-c(rownames(up_regulated),rownames(down_regulated))
length(DEG)

human<-useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
human
listAttributes(human)
attributes = c("external_gene_name","ensembl_gene_id","description","source","entrezgene_id")
human_annotation<-getBM(attributes, mart = human,uniqueRows=TRUE)
annotation<-human_annotation[human_annotation$external_gene_name%in%DEG,]


######################################
#######GO/KEGG for deg genes
#######GO/KEGG########################
BP_GO <- enrichGO(gene = annotation$entrezgene_id,
                  OrgDb= org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
BP_GO<-BP_GO@result
BP_GO
write.xlsx(BP_GO,file="GO_KEGG/BP_GO.xlsx")
BP_GO$Ontology<-"BP"

CC_GO <- enrichGO(gene = annotation$entrezgene_id,
                  OrgDb= org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
write.xlsx(CC_GO,file="GO_KEGG/CC_GO.xlsx")
CC_GO<-CC_GO@result
CC_GO$Ontology<-"CC"


MF_GO <- enrichGO(gene = annotation$entrezgene_id,
                  OrgDb= org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
MF_GO<-MF_GO@result
MF_GO$Ontology<-"MF"
write.xlsx(MF_GO,file="GO_KEGG/MF_GO.xlsx")

############################################
###Select go showing in the figures#########
############################################
head(BP_GO,30)
BP_GO<-BP_GO[c(1:8),]
head(CC_GO,30)
CC_GO<-CC_GO[c(1:8),]
CC_GO
head(MF_GO,30)
MF_GO<-MF_GO[c(1:8),]
MF_GO
###combined
GO_gp<-rbind(BP_GO,CC_GO,MF_GO)
GO_gp

# Transform the column 'Gene_number' into a numeric variable
GO_gp$Gene_number <- as.numeric(GO_gp$Count)
GO_gp$GeneRatio <- round(GO_gp$Gene_number/1079,digits = 2)
GO_gp$GeneRatio
# Replace all the "_" by a space in the column containing the GO terms
GO_gp$GO_biological_process <- chartr("_", " ", GO_gp$Description)

# Transform the column 'GO_biological_process' into factors
GO_gp$GO_biological_process<-as.factor(GO_gp$GO_biological_process)

# Transform FDR values by -log10('FDR values')
GO_gp$'|log10(FDR)|' <- -(log10(GO_gp$p.adjust))

# Change factor order
GO_gp$Group<- factor(GO_gp$Ontology,levels = c("BP","CC","MF"))
GO_gp$GO_biological_process<-factor(GO_gp$GO_biological_process,levels=GO_gp$GO_biological_process)


# Draw the plot with ggplot2
# to represent -log10(FDR) and Number of genes 
# of each GO biological process per group 
#---------------------------------------------------

GO_gp$Color<- rep(c('blue','green', 'orange'),each=8)

labs<-c("Biological Process","Cell Component","Molecular Function")
names(labs)<-c("BP","CC","MF")
GO_gp
pdf("GO_plot.pdf",width = 7,height = 6)
p<-ggplot(GO_gp, aes(x = GO_biological_process, y = GeneRatio)) +
  geom_point(data=GO_gp,aes(x=GO_biological_process, y=GeneRatio,size = Gene_number, colour = `|log10(FDR)|`,shape=Ontology), alpha=.7)+
  scale_color_gradient(low = "blue", high = "red", limits=c(0, NA))+
  coord_flip()+facet_grid(Ontology~.,scales='free',labeller=labeller(Ontology=labs))+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(color = "black"),strip.text.y = element_text(
          size = 10,face = "bold.italic"),strip.background = element_rect(colour=NA, fill="white",size=15, linetype="solid"))+
  xlab("GO enrichment")+ylab("GeneRatio")+
  labs(color="-log10(FDR)", size="Number\nof genes")

p
dev.off()




###########################
####KEGG and visualization
##########################
annotation$entrezgene_id
annotation
Kegg<- enrichKEGG(gene = annotation$entrezgene_id,organism = "human",
                  pAdjustMethod = "BH",
                  minGSSize = 1,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
dim(Kegg@result)
k<-dim(Kegg@result)[1]
k
Kegg@result
write.xlsx(Kegg@result,file="GO_KEGG/KEGG.xlsx")
Kegg_new<-Kegg@result
Kegg_new


Kegg_new<-Kegg_new[c(1:5,8:15),]
Kegg_new
Kegg_new$'|log10(FDR)|' <- -(log10(Kegg_new$p.adjust))
Kegg_new$GeneRatio <- round(Kegg_new$Count/660,digits = 2)
pdf("kegg_plot.pdf",width = 5.5,height = 5)
p<-ggplot(Kegg_new, aes(x = Description, y = GeneRatio)) +
  geom_point(data=Kegg_new,aes(x=Description, y=GeneRatio,size = Count, colour = `|log10(FDR)|`), alpha=.7)+
  scale_color_gradient(low = "blue", high = "red", limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.text = element_text(color = "black"),strip.text.y = element_text(
          size = 10,face = "bold.italic"),strip.background = element_rect(colour=NA, fill="white",size=15, linetype="solid"))+
  xlab("KEGG Pathway")+ylab("GeneRatio")+
  labs(color="-log10(FDR)", size="Number\nof genes")

p
dev.off()
