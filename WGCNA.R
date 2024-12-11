# WGCNA
rm(list = ls())
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

data <- read.csv("../all_deseq2_norm.csv",header = T,row.names = 1)
data <- data[,c(1:152)]

dir.create("../WGCNA")
setwd("./WGCNA")

wgcna.matrix <- t(data[order(apply(data,1,mad),decreasing = T)[1:5000],])
datExpr0 <-  as.data.frame(wgcna.matrix) 
datExpr0 <- log2(datExpr0+1)

# NA
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# tree
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "sample_cluster_outlier.pdf", width = 14, height = 8)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
datExpr <- datExpr0

# data traits
datTraits <- read.csv("../TEDGH_all_top5000CC/order_sample_cluster.csv",header = T,row.names = 1)
datTraits <- datTraits[rownames(datExpr),,drop = F]
all(rownames(datTraits)==rownames(datExpr))

save(datExpr, datTraits, file = "WGCNA_dataInput.RData")

# soft threshold
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf(" Analysis of network topology for various soft-thresholding powers.pdf",width = 9, height=5)
par(mfrow = c(1,2))
cex1 = 0.9 
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red") 
sft$powerEstimate
abline(h=0.90,col="red") #cutoff line
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 8 
adjacency = adjacency(datExpr, power = softPower)

# TOM
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

pdf(" Gene clustering on TOM-based dissimilarity.pdf",width = 12, height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

minModuleSize = 30    
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize) 

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf("Gene dendrogram and module colors.pdf",width = 8, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# merge similar module
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#  Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#  Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
pdf("Clustering of module eigengenes.pdf",width = 7, height=6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

pdf("Gene dendrogram and merged module colors.pdf",width = 12, height=9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepByStep.RData")

## correlation heatmap
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
colnames(datTraits)[1] <- "group"
datTraits$group <- factor(datTraits$group,levels = c("Cluster1","Cluster2"))
design <- model.matrix(~0+datTraits$group)
colnames(design) <- levels(datTraits$group)
rownames(design) <- rownames(datTraits)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) 
moduleTraitCor = cor(MEs,design, use = "p") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) 

pdf("Module-trait associations.pdf",width = 8, height=10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)  
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off()

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

design <- as.data.frame(design)
traitNames=names(design)
geneTraitSignificance = as.data.frame(cor(datExpr, design, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

y=design[,2]
GS1=as.numeric(cor(y, datExpr, use="p"))
GeneSignificance=abs(GS1)
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)
pdf(file="cluster2_GeneSignificance.pdf", width=11, height=7)
plotModuleSignificance(GeneSignificance, moduleColors)
dev.off()

# GS-MM scatter plot
trait="Cluster2"
traitColumn=match(trait,traitNames)  
for (module in modNames){
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  if (nrow(geneModuleMembership[moduleGenes,]) > 1){
    outPdf=paste(trait, "_", module,".pdf",sep="")
    pdf(file=outPdf,width=7,height=7)
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, traitColumn]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for ",trait),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    abline(v=0.8,h=0.2,col="red")
    dev.off()
  }
}

### get module genes
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

###  hub genes
geneSigFilter=0.2       
moduleSigFilter=0.8     
datMM=cbind(geneModuleMembership, geneTraitSignificance)
datMM=datMM[abs(datMM[,ncol(datMM)])>geneSigFilter,]
for(mmi in colnames(datMM)[1:(ncol(datMM)-2)]){
  dataMM2=datMM[abs(datMM[,mmi])>moduleSigFilter,]
  write.table(row.names(dataMM2), file =paste0("hubGenes_",mmi,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

# functional enrichment analysis of each module
OrgDb = "org.Hs.eg.db" 
genetype = "SYMBOL" 
table(moduleColors)
choose_module <- unique(moduleColors)

if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  tmp <- bitr(gene_module$gene,fromType = genetype, 
              toType = "ENTREZID",
              OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = OrgDb,
    ont = "BP",  
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  save(gene_module, formula_res, file="step6_module_GO_term.Rdata")
  write.csv(formula_res@compareClusterResult,
            file="step6_module_GO_term.csv")

  dotp <- dotplot(formula_res,
                  showCategory=5,
                  includeAll = TRUE,
                  label_format=100)+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(dotp,filename= "step6_module_GO_term.pdf",
         width = 16, 
         height = 15)
  
  ###run kegg analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichKEGG",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  save(gene_module, formula_res, file="step6_module_KEGG_term.Rdata")
  write.csv(formula_res@compareClusterResult,
            file="step6_module_KEGG_term.csv")

  dotp <- dotplot(formula_res,
                  showCategory=5,
                  includeAll = TRUE, 
                  label_format=100)+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  ggsave(dotp,filename= "step6_module_KEGG_term.pdf",
         width = 12, 
         height = 15)
}
