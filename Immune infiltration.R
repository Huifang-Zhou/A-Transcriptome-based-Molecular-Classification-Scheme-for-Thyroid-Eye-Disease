####  immune infiltration
rm(list = ls())
library(GSVA)
library(tibble)
library(tidyverse)
library(tidyHeatmap)
library(RColorBrewer)

expr <- read.csv("all_deseq2_norm.csv",header = T,row.names = 1)
expr <- expr[,c(1:152)]
expr <- as.matrix(log2(expr+1))

dir.create("../immuneInfiltration")
setwd("../immuneInfiltration/")

#ssGSEA
load(file = "/Users/Desktop/Public/ssGSEA28.Rdata")
gsvaP <- ssgseaParam(exprData = expr, geneSets = cellMarker,assay = NA_character_,annotation = NA_character_,
                     minSize = 1,maxSize = Inf,alpha = 0.25,normalize = T)
gsvadata <- gsva(gsvaP)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(gsvadata)

im_ssgsea <- t(ssgseaScore)
im_ssgsea <- as_tibble(im_ssgsea,rownames = "ID")
write.csv(im_ssgsea,"ssGSEA_immuneData.csv",quote = F)

ssgsea_long <- im_ssgsea %>% 
  pivot_longer(- ID,names_to = "cell_type",values_to = "Score")
head(ssgsea_long)

ssgsea_long <- left_join(ssgsea_long,sampleInf,by = "ID")
colnames(ssgsea_long)[4] <- "Cluster"

p_ssgsea <- ggplot(ssgsea_long, aes(cell_type, Score,fill = Cluster))+
  geom_boxplot(color = "black",outlier.size = 0.5,alpha = 0.8)+
  scale_fill_manual(values = mycol) + 
  theme_bw() + 
  labs(x = NULL,y = "Score") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text = element_text(color = "black",size = 8),
        legend.position = "top",
        panel.grid = element_blank())+
  stat_compare_means(aes(group = Cluster,label = ..p.signif..),
                     method = "kruskal.test",label.y = 1,size = 3)
ggsave(p_ssgsea,filename = "ssgsea_cluster_boxplot.pdf",width = 14,height = 8)

p <- ssgsea_long %>% group_by(Cluster) %>%
  heatmap(.row = cell_type,
          .column = ID,
          .value = Score,
          scale = "column",
          palette_value = circlize::colorRamp2(
            seq(-2, 2, length.out = 11), 
            rev(RColorBrewer::brewer.pal(11, "RdBu"))
          ),
          palette_grouping = list(c("#E41A1C","#377EB8")),
          show_column_names=F,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 7),
          row_title_gp = gpar(fontsize = 7),
          column_title_gp = gpar(fontsize = 7)
  )
save_pdf(p,"immune_ssgsea_heatmap.pdf")
