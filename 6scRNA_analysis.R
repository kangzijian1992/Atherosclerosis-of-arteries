library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
setwd("C:/Users/lenovo/Desktop/AS")


AS1<-read.table("RPE004_matrix.txt")
AS1<-CreateSeuratObject(AS1)
AS1$patient<-"AS1"
AS1<-RenameCells(AS1,"AS1")

AS2<-read.table("RPE005_matrix.txt")
AS2<-CreateSeuratObject(AS2)
AS2$patient<-"AS2"
AS2<-RenameCells(AS2,"AS2")

AS3<-read.table("RPE006_matrix.txt")
AS3<-CreateSeuratObject(AS3)
AS3$patient<-"AS3"
AS3<-RenameCells(AS3,"AS3")



all_data<-merge(AS1,y=c(AS2,AS3))
save(all_data,file = "All_raw_data.Rdata")

dim(all_data)
# 过滤至少检测到200个基因的细胞
selected_c <- WhichCells(all_data, expression = nFeature_RNA > 200)
# 过滤至少在3个细胞中能检测到的基因
selected_f <- rownames(all_data)[Matrix::rowSums(all_data) > 10]
length(selected_c);length(selected_f)

AS_data <- subset(all_data, features = selected_f, cells = selected_c)
dim(AS_data) 


AS_data[["percent.mt"]]<-PercentageFeatureSet(AS_data,pattern = "^mt-")
AS_data[["percent.rb"]]<-PercentageFeatureSet(AS_data,pattern = "^Rp[sl]")

VlnPlot(AS_data,features = "nFeature_RNA",group.by = "patient")
VlnPlot(AS_data,features = "nCount_RNA",group.by = "patient")

VlnPlot(AS_data,features = "percent.mt",group.by = "patient")
VlnPlot(AS_data,features = "percent.rb",group.by = "patient")

table(AS_data$percent.mt < 10  )

AS_data <- subset(AS_data, subset = percent.mt < 10 )

AS_data<-AS_harmony
DefaultAssay(AS_data) <- "RNA"
AS_data <- NormalizeData(AS_data)

AS_data <- FindVariableFeatures(AS_data, selection.metASd = "vst", nfeatures =2000)

AS_data <- ScaleData(AS_data,vars.to.regress = "percent.mt")

AS_data <- RunPCA(AS_data,features=VariableFeatures(object = AS_data),npcs = 50)

library(harmony)

AS_harmony<- AS_data %>% RunHarmony("patient",plot_convergence=T)

harmony_embeddings <- Embeddings(AS_harmony,"harmony")


AS_harmony<- AS_harmony%>% RunUMAP(reduction="harmony",dims=1:30) %>%
  FindNeighbors(reduction="harmony",dims=1:30) %>%
  FindClusters(resolution=2) %>%
  identity()

p1<-DimPlot(object = AS_harmony, reduction = "umap", group.by = "patient")
p2<-DimPlot(object = AS_harmony, reduction = "umap",label = T)

## Annotation

AS_genelist<-c("CD3D","CD79A","CD79B","MS4A1","IGHM","IGHD","CD27","CD38","CD80","CD86",
               "SLAMF7","SDC1","PRDM1","MZB1","IGHG1","IGHG2","IGHG3","NOTCH2","JAM3",
               "CCR7","IGHA1","TNFSF17","CXCR4","CXCR5","CR2","BANK1","ITGA4","ITGAX")

AS_genelist<-c("LY6A","LY6C1","VCAM1","LYZ","PDGFRB","RSG5","CD14","FCGR3A","CD163","COL1A2",
               "CD1C","PECAM1","SLC38A","VWF","COL4A1","IL1B")

AS_genelist<-c("CD79A","CD3D","MPO","TPSB2","DCN","RGS5","ABCC9","ACTA2","SOX17")
genelist<-intersect(rownames(AS_harmony),AS_genelist)

mytheme2 <- theme(plot.title = element_text(face = "bold.italic",
                                            size = "30", color = "brown"),
                  axis.title = element_text(face = "bold.italic",
                                            size = "20",color = "blue"),
                  axis.text.x = element_text(face = "bold",
                                             size = 20, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(face = "bold",size = 10),
                  panel.background = element_rect(fill = "white", 
                                                  color = "black"),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 20,
                                              face = "bold"),
                  panel.grid.minor.x = element_blank())

for(i in 1:length(genelist)){
  name <- paste(genelist[i],'plot.png',sep='_')
  png(name,width=800,height=800)
  gg <- FeaturePlot(AS_harmony,features =genelist[i])+
    # NoLegend()+
    mytheme2
  
  print(gg)
  dev.off()
}

AS_harmony<-RenameIdents(AS_harmony,"3"="T_cell",
                          "18"="B_cell",
                          "14"="Myeloid",
                          "9"="Myeloid",
                          "17"="Myeloid",
                          "2"="Myeloid",
                          "8"="Myeloid",
                          "13"="Mast_cell",
                          "15"="SMC",
                          "4"="SMC",
                          "16"="SMC",
                          "5"="SMC",
                          "0"="SMC",
                          "6"="SMC",
                          "1"="SMC",
                          "10"="SMC",
                          "11"="Endo1",
                          "7"="Endo1",
                          "12"="Endo2"
                          
)
table(AS_harmony$patient,AS_harmony@active.ident)
DimPlot(AS_harmony)
save(AS_harmony,file = "AS_harmony.Rdata")

## UMAP plot

gg<-DimPlot(AS_harmony,label = F)
gg
ggsave(gg,filename = "AS_UMAP.png",width = 6,height = 4)

## Dot plot

marker<-FindAllMarkers(AS_harmony)
ribo.genes <- unique(grep(pattern = "^RP[SL]", x = marker$gene, value = TRUE))
mito.genes <- unique(grep(pattern = "^MT-", x = marker$gene, value = TRUE))
x<-ifelse(marker$gene %in% c(mito.genes,ribo.genes),FALSE,TRUE)
marker <- marker[x,]
write.csv(marker,file ="AS_marker.csv" )
top5 <- marker %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
top5_gene <- as.character(top5$gene)
top5_gene<-unique(top5_gene)

gg<-DotPlot(AS_harmony, features =top5_gene , cols = "RdYlBu") + RotatedAxis()
gg
ggsave(gg,filename = "AS_dot.png",width = 16,height = 6)

## Dotplot pathological gene


library(Seurat)
library(ggplot2)
DimPlot(AS_harmony)

gene<-read.csv("keygenes.csv")
gene<-gene$x
gene_all<-as.character(gene)

gene1<-c("MDM2","GNAI2","NCSTN","PARVB","XPO6","ADAP2","PLEKHO2","CD68",
         "NCKAP1L","SLC7A8","SLC3A2","MAF","LY6E","ADM","FABP4","TDO2")


gene2<-gene_all[-match(gene1,gene_all)]

genelist<-list("unstable gene"=gene1,"stable gene"=gene2)


cluster<-c("SMC","SEM","FC","Fibro1","Fibro2","Endo_Artery","Endo_Capillaries","Monocyte",
           "M1_Macro","M2_Macro","DC","CD4 T","CD8 T","NK","B","Mast")
AS_harmony$cluster<-AS2_merge_harmony@active.ident
AS_harmony$cluster<-factor(AS_harmony$cluster,levels = rev(cluster))
gg<-DotPlot(AS_harmony,features = genelist,cols = "RdYlBu",group.by = "cluster",dot.scale = 4.5)+theme_classic()+
  RotatedAxis()+theme(axis.text.x.bottom = element_text(colour = "black"),axis.text.y.left = element_text(colour = "black"),
  )
ggsave(gg,filename = "stable_unstable_gene.png",width = 11,height = 4)



## Stability score



AS_harmony<-AddModuleScore(AS_harmony,features = list(gene1,gene2),name = c("unstable","stable"))

head(AS_harmony@meta.data)

table(AS_harmony@active.ident)


cols<-c("SMC"="#FF68A1","SEM"="#FF61CC","FC"="#00BE67",
        "Fibro1"="#00C19A","Fibro2"="#00BFC4","CD4 T"="#E68613",
        "CD8 T"="#CD9600","NK"="#ED68ED","B"="#0CB702","Mast"="#8494FF","M1_Macro"="#00B8E7","M2_Macro"="#00A9FF","Monocyte"="#C77CFF",
        "DC"="#ABA300","Endo_Artery"="#7CAE00","Endo_Capillaries"="#0CB702")

gg1<-VlnPlot(AS_harmony,features ="unstale1",pt.size = 0)+scale_x_discrete(limits=c("SMC","SEM","FC",
                                                                                          "Fibro1","Fibro2","CD4 T",
                                                                                          "CD8 T","NK","B","Mast","M1_Macro","M2_Macro","Monocyte",
                                                                                          "DC","Endo_Artery","Endo_Capillaries"))+NoLegend()+
  scale_fill_manual(values = cols)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",size=0.65,
               width = 0.25,position = position_dodge( .9))+
  stat_summary(fun.y = mean,color = "white",geom = "point",size = 1)

gg2<-VlnPlot(AS_harmony,features ="stable2",pt.size = 0)+scale_x_discrete(limits=c("SMC","SEM","FC",
                                                                                   "Fibro1","Fibro2","CD4 T",
                                                                                   "CD8 T","NK","B","Mast","M1_Macro","M2_Macro","Monocyte",
                                                                                   "DC","Endo_Artery","Endo_Capillaries"))+NoLegend()+
  scale_fill_manual(values = cols)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",size=0.65,
               width = 0.25,position = position_dodge( .9))+
  stat_summary(fun.y = mean,color = "white",geom = "point",size = 1)

ggsave(gg1,filename = "unsatble_score.png",width = 6,height = 4.5)
ggsave(gg2,filename = "satble_score.png",width = 6,height = 4.5)