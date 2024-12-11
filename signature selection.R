# signature selection
rm(list = ls())
library(glmnet)
library(ggpubr)

data <- read.csv("../all_deseq2_norm.csv",header = T,row.names = 1)
data <- data[,c(1:152)]
d <- log2(data+1)

greenyellow_hub <- read.table("hubGenes_MMgreenyellow.txt")
greenyellow_hub <- greenyellow_hub$V1
green_hub <- read.table("hubGenes_MMgreen.txt")
green_hub <- green_hub$V1

deg_intersect <- c(greenyellow_hub,green_hub)

d <- d[deg_intersect,]

load("order_sample_cluster_matrix_anno_forHeatmap.Rdata")
d <- as.data.frame(t(d))
all(rownames(d)==rownames(annCol5))
d <- d[rownames(annCol5),]
d$cluster <- annCol5$Cluster

set.seed(1234)

x <- as.matrix(d[,1:(ncol(d)-1)])
y <- as.matrix(d$cluster)

alpha1_fit <- glmnet(x, y, alpha = 1, family = "binomial", nlambda = 100)

pdf("lasso.pdf",width = 6,height = 6)
plot(alpha1_fit, xvar = "lambda", label = TRUE)
alpha1.fit.cv <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")
plot(alpha1.fit.cv)
dev.off()

print(alpha1.fit.cv)
coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)

feature_all <- as.data.frame(as.matrix(coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)))
colnames(feature_all) <- "coff"
feature_opt <-  feature_all %>% filter(abs(coff) > 0)
write.table(rownames(feature_opt),"lassogene.txt")

# key gene boxplot
data <- log2(data+1)

gene <- read.table("lassogene.txt")
gene <- gene$x
gene <- gene[-1]

data <- data[rownames(data) %in% gene,]
load("order_sample_cluster_matrix_anno_forHeatmap.Rdata")
data <- as.data.frame(t(data))
data <- data[rownames(annCol5),]
data$cluster <- annCol5$Cluster

rt <- data
mycol = c("Cluster1" ="#E41A1C","Cluster2" = "#377EB8")
dir.create("./lassogene_Boxplot")

for(gene in colnames(rt[,1:ncol(rt)-1])){
  data=rt[,c(gene,"cluster")]
  group=levels(factor(data$cluster))
  data$cluster=factor(data$cluster, levels=group)
  comp=combn(group,2) 
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  boxplot=ggboxplot(data, x = "cluster", y= gene,fill = "cluster",
                    xlab="Cluster",
                    ylab=paste0(gene," Expression\nLog2 (Normalize+1)"),
                    legend.title= "Cluster",
                    alpha = 0.8)+ 
    geom_jitter(position = position_jitter(0.2),size=0.5)+
    stat_compare_means(comparisons = my_comparisons)+
    scale_fill_manual(values = mycol)
  
  pdf(file=paste0("./lassogene_Boxplot/", gene, ".pdf"), width=4, height=5)
  print(boxplot)
  dev.off()
}
