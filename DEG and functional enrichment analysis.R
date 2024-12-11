# Differential expression and functional enrichment analysis
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

data <- read.csv("/Users/Desktop/TAO/blood_bulk_RNAseq/ERCC_normalized_count_matrix.csv",header = T,row.names = 1)

# data preprocessing
keep <- rowSums(data>0) >= floor(0.75*ncol(data))
table(keep)
filter_count <- data[keep,]

# sample group
sampleinf <- data.frame(sample = colnames(filter_count),condition = c(rep("TED",152),rep("GH",20)))
sampleinf$condition <- factor(sampleinf$condition,levels = c("TED","GH"))

# DEG
data <- filter_count
dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, 
                                      colData = sampleinf, design = ~ condition)
vsd <- vst(dds,blind = FALSE)
boxplot(assay(vsd),las = 2)
DESeq2::plotPCA(vsd)
# PCA
library(factoextra)
library(FactoMineR)
dat <- as.data.frame(t(assay(vsd)))
dat_pca <- PCA(dat, graph = FALSE)
p <- fviz_pca_ind(dat_pca,
                   geom.ind = "point", 
                   col.ind = sampleinf$condition, 
                   palette = c("#00AFBB", "#E7B800","#EF6F6AFF"), 
                   addEllipses = T,  
                   legend.title = "Group")+
  ggtitle("PCA")+theme(plot.title = element_text(size = 12,hjust = 0.5))
ggsave(p,filename = "PCA_Group.pdf")

# normalization and wald test
dds <- DESeq2::DESeq(dds,test="Wald")
res <- DESeq2::results(dds, contrast = c("condition", "TED", "GH")) %>% as.data.frame() 

# 设置过滤p值和log2FC，寻找DEG
Diff_FC <- 1.5
pvalue <- 0.05
res$Class <- ifelse(res$log2FoldChange > log2(Diff_FC) & res$pvalue < pvalue,"Up",ifelse(res$log2FoldChange < -log2(Diff_FC) & res$pvalue < pvalue,"Down","None" ))

# normalized matrix 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
write.csv(normalized_counts,"all_deseq2_norm.csv",quote = F)

DEseqExp <- res

# remove NA
DEseqExp <- DEseqExp[!is.na(DEseqExp$log2FoldChange),]
DEseqExp$Gene <- row.names(DEseqExp)

# save final results
save(dds,file = "all_dds.Rdata")
dir.create("./TED_VS_GH")
write.csv(DEseqExp, "./TED_VS_GH/TED_VS_GH_FC_1.5_pvalue_0.05_DEG.csv",quote = F,row.names = T)
#DEG over-----------

# volcano plot
setwd("./TED_VS_GH/")
DEseqExp <- read.csv("TED_VS_GH_FC_1.5_pvalue_0.05_DEG.csv",header = T,row.names = 1)
sub <- DEseqExp[!is.na(DEseqExp$Class),]
sub$pvalue_L <- -log10(sub$pvalue+10^-300)

sub_Labels1 <- sub %>% dplyr::filter(Class %in% c("Up","Down") ) %>% dplyr::group_by(Class)  %>% 
  dplyr::top_n(10,pvalue_L)

label <- c("TSHR","IGF1R","ERBB4","KDR", "EGFR","FAT1", "MST1R", "FGFR2", "NTRK2", "AKAP13",
           "NTRK3","EPHA7", "PTPRT", "EPHA3", "NLRP1","CD44", "FCGR3B","IL1R2", "IL23R","IL17C")
sub_Labels2 <- sub[sub$Gene %in% label,]
sub_Labels <- rbind(sub_Labels1,sub_Labels2)

volcano_p <- ggplot(sub,aes(x=log2FoldChange ,y= pvalue_L))+
  geom_point(aes(color=Class),size=0.8)+
  labs(x="log2(Foldchange)",y="-log10(p-value)", title = "TED_VS_GH")+
  scale_color_manual(limits=c("Up","Down","None"),values=c("#e41a1c","#3182bd","lightgray"),#values=c(NA,NA,NA),#
                     labels=c(paste("Up (n=",table(sub$Class)[3],")",sep=""),
                              paste("Down (n=",table(sub$Class)[1],")",sep=""),
                              paste("None (n=",table(sub$Class)[2],")",sep="")),name="")+
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)),linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme(panel.background = element_rect(fill=NA,color="black"),axis.title= element_text(color="black",size=12),panel.grid.major = element_blank(),
        axis.text= element_text(color="black",size=10),legend.position = "bottom",
        legend.background = element_blank(),legend.key = element_blank())+
  ggrepel::geom_label_repel(data = sub[sub$Gene %in% sub_Labels$Gene, ],
                            aes(label = Gene),
                            force = 1,
                            size = 3,               
                            fill = "white",         
                            color = "black",       
                            label.size = 0.25,     
                            box.padding = 0.35,    
                            label.padding = 0.2)   
ggsave("DEG_FC1.5_volcano_plot_label.pdf",volcano_p,width = 8,height = 8)

# get DEGs
up <- rownames(DEseqExp)[DEseqExp$Class == "Up"]
down <- rownames(DEseqExp)[DEseqExp$Class == "Down"]
diff <- c(up,down)

# heatmap
normalized_counts <- read.csv("../all_deseq2_norm.csv",header = T,row.names = 1)
hmExp <- normalized_counts[rownames(normalized_counts)%in% label,]
hmExp <- log2(hmExp+1)
hmExp <- hmExp[label,]

# annotation
Type <- data.frame(row.names = colnames(hmExp),group = c(rep("TED",152),rep("GH",20)))
annorow <- data.frame(row.names = label,Category = c(rep("Autoimmune antigens",2),rep("Growth factor related genes",6),
                                                     rep("Nervous system related genes",6),rep("Immune response related genes",6)))

annColors <- list(group = c("TED" ="#E64B35B2","GH" = "#00A087B2" ),
                  Category =c("Autoimmune antigens" = "#FDB462","Growth factor related genes" = "#FFEDA0",
                              "Nervous system related genes"="#80B1D3","Immune response related genes"="#8DD3C7"))

pdf(file="/Users/qianweijin/Desktop/TAO/blood_bulk_RNAseq/20240916TED_GH/TED_VS_GH/heatmap_label.pdf",height=5,width=6)
pheatmap(hmExp, 
         scale = "row",
         annotation_col=Type, 
         annotation_row = annorow,
         annotation_colors = annColors,
         color = colorRampPalette(colors = c("#2166ac","#f7fbff","#b2182b"))(50),
         cluster_cols =F,
         cluster_rows = F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=8,
         fontsize_col=10)
dev.off()


# functional enrichment analysis
# GO
go_up <- enrichGO(up_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff = 0.05,qvalueCutoff = 0.1,
                  pAdjustMethod = "BH", 
                  minGSSize = 10,maxGSSize = 500,readable = T)
go_res_up <- go_up@result
write.csv(go_res_up,"GO_up.csv")

go_down <- enrichGO(down_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pvalueCutoff = 0.05,qvalueCutoff = 0.1,
                    pAdjustMethod = "BH", 
                    minGSSize = 10,maxGSSize = 500,readable = T)
go_res_down <- go_down@result
write.csv(go_res_down,"GO_down.csv")

# KEGG
up_entrez <- bitr(up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
KEGG_up <- enrichKEGG(up_entrez$ENTREZID,organism = "hsa", 
                      pvalueCutoff = 0.05,qvalueCutoff = 0.1,
                      pAdjustMethod = "BH", 
                      minGSSize = 10,maxGSSize = 500)
KEGG_up <- setReadable(KEGG_up,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
KEGG_res_up <- KEGG_up@result
write.csv(KEGG_res_up,"KEGG_up.csv") 

down_entrez <- bitr(down,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
KEGG_down <- enrichKEGG(down_entrez$ENTREZID,organism = "hsa",  
                        pvalueCutoff = 0.05,qvalueCutoff = 0.1,
                        pAdjustMethod = "BH", 
                        minGSSize = 10,maxGSSize = 500)
KEGG_down <- setReadable(KEGG_down,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
KEGG_res_down <- KEGG_down@result
write.csv(KEGG_res_down,"KEGG_down.csv") 

# GSEA
h <- read.gmt("/Users/Desktop/Public/pathway_gmt/h.all.v2023.2.Hs.symbols.gmt")
h$term <- gsub("HALLMARK_","",h$term)
genelist <-  DEseqExp$log2FoldChange
names(genelist) <-  rownames(DEseqExp)
genelist = sort(genelist,decreasing = T) 
h_ges <- GSEA(genelist,
              TERM2GENE = h,
              minGSSize = 10,
              maxGSSize = 500,
              pvalueCutoff = 0.5,
              pAdjustMethod = "BH",
              verbose = F,
              eps = 0)
h_res <- h_ges@result
write.csv(h_res,"GSEA_Hallmark_result.csv",quote = F)

cp_kegg <- read.gmt("/Users/Desktop/Public/pathway_gmt/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
cp_reactome <- read.gmt("/Users/Desktop/Public/pathway_gmt/c2.cp.reactome.v2023.2.Hs.symbols.gmt")
c5_gobp <- read.gmt("/Users/Desktop/Public/pathway_gmt/c5.go.bp.v2023.2.Hs.symbols.gmt")

kegg_ges <- GSEA(genelist,
                 TERM2GENE = cp_kegg,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = F,
                 eps = 0)
kegg_res <-kegg_ges@result
write.csv(kegg_res,"GSEA_cpKEGG_result.csv",quote = F)

react_ges <- GSEA(genelist,
                  TERM2GENE = cp_reactome,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  verbose = F,
                  eps = 0)
react_res <-react_ges@result
write.csv(react_res,"GSEA_cpREACTOME_result.csv",quote = F)

gobp_ges <- GSEA(genelist,
                 TERM2GENE = c5_gobp,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose = F,
                 eps = 0)
gobp_res <-gobp_ges@result
write.csv(gobp_res,"GSEA_c5GOBP_result.csv",quote = F)

save(h_ges,kegg_ges,react_ges,gobp_ges,file = "GSEA_ges.Rdata")
