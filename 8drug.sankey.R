drug_genes<-read.delim("interactions.tsv",header = T,sep = "\t")
head(drug_genes)
summary(drug_genes$interaction_group_score)
drug_genes<-drug_genes[!drug_genes$gene_name=="NA",]
drug_genes<-drug_genes[!drug_genes$gene_name=="",]
drug_genes
#drug_genes<-drug_genes[drug_genes$interaction_group_score>1,]



key_genes<-read.xlsx("keygenes.xlsx",sheetIndex = 1)
key_genes$x
#key_genes<-c("CNN1","IGFBP6","TMEM43","TMEM47")
drug_genes<-drug_genes[drug_genes$gene_name%in%key_genes$x,]
drug_genes
drug_genes$interaction_group_score

drug_genes<-drug_genes[!is.na(drug_genes$interaction_group_score),] 
drug_genes<-drug_genes[drug_genes$interaction_group_score,]
write.xlsx(drug_genes,file="drug_genes.xlsx")

graph
colnames(graph)[1:2]<-c("immune_cells","gene_name")
imuu.cell<-graph[,1:2]


drug.genes<-drug_genes[,c("gene_name","drug_claim_name")]


library(reshape)
library(ggalluvial)
#install.packages("ggalluvial")

interaction <- merge(imuu.cell, drug.genes, by = 'gene_name')
interaction
interaction$link <- 1
interaction <- reshape::melt(interaction, id = 'link')
interaction
variable <- summary(interaction$variable)



interaction$flow <- rep(1:variable[1], length(variable))
head(interaction)


interaction
mycol <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462',
           
           '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
           
           '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5',
           
           '#6181BD', '#F34800', '#64A10E', '#FF00FF', '#c7475b', '#049a0b', '#BEAED4',
           
           '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#4253ff', '#ff4308', '#D8D155',
           
           '#64495D', '#7CC767',"orange","green","skyblue","gray")

interaction
p <- ggplot(interaction, aes(x = variable, y = link,
                       
                       stratum = value, alluvium = flow, fill = value)) +
  
  geom_stratum() + #冲击图中的堆叠柱形图
  geom_flow(aes.flow = 'forward') + #冲击图连线绘制
  
  scale_fill_manual(values = mycol) + #颜色赋值
  
  geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) + #添加 lncRNA、miRNA 和 mRNA 标签
  
  scale_x_discrete(limits = c('drug_claim_name', 'gene_name', 'immune_cells')) + #定义 lncRNA、miRNA 和 mRNA 列的展示顺序
  
  labs(x = '', y = '') + #去除 x 轴和 y 轴标题
  
  theme(legend.position = 'none', panel.background = element_blank(),
        
        line = element_blank(), axis.text.y = element_blank()) #去除背景和图例

print(p)
pdf("/Users/yaoyuelin/Desktop/drug.pdf")
print(p)
dev.off()

