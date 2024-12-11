### cMap
rm(list = ls())
library(data.table)
library(circlize)
library(ComplexHeatmap)
library(pheatmap)

setwd("../drug")
setwd("../cluster2")
res <- fread("cluster2 green.txt",,data.table = F)
head(res)

res_sorted <- res[order(res$Score, decreasing = TRUE), ]
res_sorted <- dplyr::filter(res_sorted,res_sorted$Description != "-666")
res_tail50 <- tail(res_sorted[,c(2,5,6)],50)
write.csv(res_tail50,"res_tail50.csv")

res_tail50 <-res_tail50[,c("Name","Score","Description")]
rownames( res_tail50)<-NULL
mydata<- res_tail50[,-3]%>%
  column_to_rownames('Name')

input_raw <- res[res$Name%in%rownames(mydata),c(6,5)]
input_raw <- input_raw[!duplicated(input_raw$Name),]
input <- reshape2::dcast(input_raw, Description ~ Name)
rownames(input) <- input$Description 
input <- input[, -1] 
input <- as.matrix(input)

input[!is.na(input)] <- "inhibitor"
input[is.na(input)] <- ""
input <- input[, order(colnames(input))]

pp <- as.data.frame(input)
pp <- data.frame(Name=colnames(pp),Type=1:50)
dd <- merge(input_raw,pp,by.x = 2,by.y=1)
dd$Type[dd$Description%like%'inhibitor'] <- "Inhibitor"
dd$Type[dd$Description%like%'agonist'] <- "Agonist"
dd$Type[dd$Description%like%'channel blocker'] <- "Channel blocker"
dd$Type[dd$Description%like%'exchange inhibitor'] <- "Exchange inhibitor"
dd$Type[dd$Description%like%'analog'] <- "Analog"
dd$Type[dd$Description%like%'antagonist'] <- "Antagonist"
dd$Type[dd$Type %in% c(0:50)] <- "Others"
dd<- dd[order(dd$Type,decreasing = F),]
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "lightblue", col = "grey")),
  # dots
  inhibitor = function(x, y, w, h) 
    grid.points(x, y, pch = 16, size = unit(0.8, "char"))
  
)

ha_rowdata <- rowSums(apply(input, 2, function(x) x=="inhibitor") + 0) %>% as.numeric() 
table(dd$Type)
dd$Type <- factor(dd$Type,levels = c("Agonist","Antagonist",
                                     "Inhibitor","Others"))
top_ha <- HeatmapAnnotation(Type = dd[,3], 
                            col = list(Type = c("Agonist" = "#FDBF6F",#1666A5
                                                "Antagonist" = "#33A02C",
                                                "Channel blocker" = "#1F78B4",#B5CEE0
                                                "Others" = "#6A3D9A",
                                                "Inhibitor"="#E31A1C")),#F7E897FF
                            show_legend = T,
                            annotation_height = unit(30, "mm"),
                            #annotation_legend_param = list(Type = list(title = "Type")),
                            annotation_name_side = "left",
                            annotation_name_rot = 90)
right_ha <- rowAnnotation(count = anno_barplot(ha_rowdata, axis = F, border = F, 
                                               gp = gpar(fill = "lightblue"),
                                               bar_width = 1, width = unit(2, "cm")),
                          annotation_name_side = "top",
                          annotation_name_rot = 0)


column_order <- intersect(dd$Name,colnames(input))
input <- input[,column_order]
oncoPrint(input, alter_fun = alter_fun, 
          show_column_names = TRUE, 
          column_names_side = "top",
          column_order = column_order, 
          top_annotation = top_ha,
          right_annotation = right_ha,
          show_pct = F, 
          show_heatmap_legend = F)
decorate_annotation("Type", {
  grid.text("Mechanism of Action", unit(1, "npc") + unit(3, "mm"), just = "left")})

graph2pdf(file=paste0('oncoPrint_cluster2.pdf'),width = 20,height = 10)
dev.off()
