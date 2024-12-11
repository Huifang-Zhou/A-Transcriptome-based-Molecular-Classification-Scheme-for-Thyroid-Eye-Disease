# Consensus Cluster
rm(list = ls())
library(ggvenn)
library(ConsensusClusterPlus)
library(limma)
library(GSEABase)
library(ggplot2)
library(WGCNA)
library(ggsci)
library(RColorBrewer)
library(pheatmap)

data <- read.csv("all_deseq2_norm.csv",header = T,row.names = 1)
d <- data[,1:152]
d <- log2(d+1)

mads=apply(d,1,mad)
d=d[rev(order(mads))[1:5000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

tao_topgene <- rownames(d)
dir.create("./TEDGH_all_top5000CC/")
write.table(tao_topgene,"./TEDGH_all_top5000CC/TED_TOP5000_gene.txt")

deg <- read.csv("./TED_VS_GH/TED_VS_GH_FC_1.5_pvalue_0.05_DEG.csv",header = T,row.names = 1)
deg <- rownames(deg[deg$Class %in% c("Up","Down"),])

# gene intersection
deg_intersect <- Reduce(intersect, list(deg,tao_topgene))
write.table(deg_intersect,"./TEDGH_all_top5000CC/tedgh_all_top5000_intersect_gene_for_cc.txt")

# venn plot
mycol <- c("#E64B3566", "#4DBBD566" )
deg <- list("TED VS GH" = deg,"TED TOP 5000 GENE" = tao_topgene)
p_deg <- ggvenn(deg,show_percentage = F,
                stroke_color = "white",
                fill_color = mycol)
ggsave(p_deg,filename = "./TEDGH_all_top5000CC/new_intersect_venn_plot.pdf")


# Consensusclusterplus-----------
data <- as.matrix(d[deg_intersect,])
setwd("./TEDGH_all_top5000CC/")
workDir = getwd()
maxK=9
results = ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=50, 
                               pItem=0.8, 
                               pFeature=1,
                               title=workDir,
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="pdf")

Kvec = 2:9
x1 = 0.1; x2 = 0.9       
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")  
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])        
  PAC[i-1] = Fn(x2) - Fn(x1)
} 
optK = Kvec[which.min(PAC)] 
optK

data <- d[deg_intersect,]
datExpr0=t(data)

##############
clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]
cluster = paste0("Cluster",cluster)
cluster <- data.frame(row.names = colnames(data),cluster = cluster)
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)

# pca analysis
data= t(data)
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)

cluster=read.table("cluster.txt",sep="\t",header=F)               
group=paste0("",as.vector(cluster[,2]))
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)

pdf(file="PCA.pdf",height=5,width=6.5)     
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

###########
colnames(cluster)[2] = "Cluster"
# annotation
annCol <- read.csv("../sample_anno.csv",header = T,row.names = 1)
all(rownames(annCol) == cluster$V1)
annCol = cbind(annCol, cluster)
annCol = annCol[,-which(colnames(annCol)%in%"V1")]
plotdata <- as.data.frame(t(data))
dim(plotdata)
dim(annCol)

# arrange
plotdata2 = cbind(t(plotdata),annCol)
dim(plotdata2)
plotdata3 = plotdata2[order(plotdata2$Cluster),]
sample_cluster_inf <- plotdata3[,ncol(plotdata3)] %>% as.data.frame()
rownames(sample_cluster_inf) <- rownames(plotdata3)
colnames(sample_cluster_inf) <- "sample_cluster"
write.csv(sample_cluster_inf,file="order_sample_cluster.csv",quote=F,row.names = T)  

sample_cluster_inf$sample <- rownames(sample_cluster_inf)
num <- as.data.frame(strsplit(sample_cluster_inf$sample,"_"))
num <- as.data.frame(t(num))
num <- num$V3
sample_cluster_inf$sample <- as.numeric(num)

order_sample_cluster <- arrange(sample_cluster_inf,sample)
write.csv(order_sample_cluster,file="order_sample.csv",quote=F,row.names = T)  

a = nrow(plotdata)+1
b = ncol(plotdata3)

plotdata4 = plotdata3[,-(a:b)]

plotdata5 = t(plotdata4)

annCol5 =  plotdata3[,-(1:a-1)]

save(plotdata5,annCol5,file = "order_sample_cluster_matrix_anno_forHeatmap.Rdata")

annColors <- list(Cluster = c("Cluster1" = "#E41A1C","Cluster2" =  "#377EB8"),
                  gender = c("male" = "#8DA0CB", "female" = "#E78AC3"),
                  activity = c("active" = "#FC8D62", "inactive" = "#66C2A5"),
                  state = c("mild" = "#FFEDA0", "moderate to severe" = "#FEB24C", "very severe" = "#F03B20"))

# heatmap
set.seed(123)
pdf(file = "cluster_heatmap_order.pdf",width = 20, height = 12)
pheatmap(plotdata5,
         scale = "row",
         border_color = NA,
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = T,
         show_rownames = T,
         annotation_col = annCol5,
         annotation_colors = annColors,
         gaps_col =  cumsum(table(annCol5$Cluster))[1:2], 
         fontsize = 12,
         fontsize_col = 4,
         fontsize_row = 1,
         color = colorRampPalette((c("#2166AC","white","#D53E4F")))(64)) 
dev.off()

plotdata6 <- t(scale(t(plotdata5)))
max(plotdata6)
min(plotdata6)
plotdata6[plotdata6 > 4] <- 4
plotdata6[plotdata6< -4] <- -4

pdf(file = "cluster_heatmap_order.pdf",width = 20, height = 12)
pheatmap(plotdata6,
         border_color = NA,
         cluster_cols = F,
         cluster_rows = T,
         show_colnames = T,
         show_rownames = T,
         annotation_col = annCol5,
         annotation_colors = annColors,
         gaps_col =  cumsum(table(annCol5$Cluster))[1:2],  
         fontsize = 12,
         fontsize_col = 4,
         fontsize_row = 0.5,
         color = colorRampPalette((c("#2166AC","white","#D53E4F")))(64)) 
dev.off()

