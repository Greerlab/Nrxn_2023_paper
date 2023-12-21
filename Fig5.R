library(Seurat)
library(harmony)
library(patchwork)
library(ggplot2)
library(viridis)
library(pheatmap)
library(reshape2)
MTC = readRDS("data/MTC.rds")
vizgen.obj = readRDS("data/Sample_0601_filtered_50_log.rds")

genes_sjp = c("Nrxn1","Nrxn2","Nrxn3",
              "Nxph1","Nxph2","Nxph3","Nxph4",
              "Nlgn1","Nlgn2","Nlgn3",
              "Lrrtm1","Lrrtm2","Lrrtm3","Lrrtm4",
              "Cbln1","Cbln2","Cbln3","Cbln4",
              "Adgrl1","Adgrl2","Adgrl3","Adgrl4",
              "C1ql1","C1ql2","C1ql3", "C1ql4",
              "Clstn1","Clstn2","Clstn3",
              "Dag1","Car10","Car11",
              "Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6",
              "Tafa1","Tafa2","Tafa3","Tafa4","Tafa5")
genes_sjp1 = grep("Nxph|Nlgn|Lrtm|Cbln|Lphn|C1ql|Clstn|Gabra|Dag",genes_sjp, value = T)

# 5A
em = as.matrix(GetAssayData(MTC, slot = "data")[genes_sjp1,])
em = reshape2::melt(em)
ggplot(em, aes(Var1, value, fill = Var1))+geom_violin(scale = "width")+theme_bw()

# 5B
DimPlot(MTC, group.by = "MTC_sub_group", label = T)+NoLegend()

# 5C
Idents(MTC)

MTC = ScaleData(MTC, features = rownames(MTC))
MTC_ave = AverageExpression(MTC, return.seurat = T)
mtx = GetAssayData(MTC_ave, slot = "scale.data")[genes_sjp1,]
mtx[mtx>=2] = 2
mtx[mtx<= (-2)] = (-2)
p = DoHeatmap(MTC_ave, features = genes_sjp, draw.lines = F)
myColor <- colorRampPalette(c("blue", "white", "Red"))(50)
myBreaks <- c(seq(min(mtx), 0, length.out = 26), seq(0.1, max(mtx), , length.out = 24)) 
pheatmap(mtx)

# 5E
vizgen.obj = readRDS("~/Dropbox (UMass Medical School)/OB_MERFISH/220601/Sample_0601_filtered_50_log.rds")
Idents(MTC) = MTC$MTC_sub_group

vizgen.obj = subset(vizgen.obj, idents = "MCL")
anchors <- FindTransferAnchors(reference = MTC, query = vizgen.obj, normalization.method = "LogNormalize",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = MTC$MTC_sub_group, prediction.assay = TRUE,
                                  weight.reduction = vizgen.obj[["pca"]], dims = 1:30)
vizgen.obj[["predictions"]] <- predictions.assay
DefaultAssay(vizgen.obj) <- "predictions"
clusters = unique(MTC$MTC_sub_group)

ImageFeaturePlot(vizgen.obj,features = "ET1",fov = "OB1")
plots = c()
for (i in 1:length(clusters)) {
  plots[[i]] = ImageFeaturePlot(vizgen.obj,features = clusters[i],fov = "OB1")+NoLegend()+ggtitle(paste0("cluster ",clusters[i]," projection"))&scale_fill_viridis(direction = 1,option = "C")
}
plots[[length(clusters)+1]] = DimPlot(MTC, label = T, label.size = 4)+NoLegend()
p = wrap_plots(plots, nrow = 2)
p

tmp = as.data.frame(vizgen.obj@assays$predictions@data)
predicted_cluster = rownames(tmp)[apply(tmp, 2, which.max)]
predicted_cluster[apply(tmp, 2, max)<0.4]= "unclassified"
vizgen.obj$predicted_cluster = factor(predicted_cluster)
p = ImageDimPlot(vizgen.obj,group.by = "predicted_cluster" ,fov = "OB1", axes = TRUE, size = 0.5, cols = "polychrome",dark.background = T, na.value = "black")
p