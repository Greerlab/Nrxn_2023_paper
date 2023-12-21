library(Seurat)
library(org.Mm.eg.db)
library(GO.db)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(pheatmap)

# 1A-C
OSN = readRDS("data/OR_OSN.rds")
OSN = AverageExpression(OSN, return.seurat = T)
OSN = FindVariableFeatures(OSN)
k = keys(org.Mm.eg.db, keytype="GO")
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("ONTOLOGY"), keytype="GO")
df = df[df$ONTOLOGY=="BP",]
BPGO = df$GO
k = unique(BPGO)
df = AnnotationDbi::select(org.Mm.eg.db, keys=k, columns=c("SYMBOL"), keytype="GO")
df = df[,c("GO","SYMBOL")]
df = df[df$SYMBOL%in%rownames(OSN),]%>% distinct()

x = df %>%
  group_by(GO) %>%
  summarise(count = n_distinct(SYMBOL))
x = x[x$count>=150,]
ref = x$GO

goterms <- cbind.data.frame("GO" = names(Term(GOTERM)), "GO_name" = Term(GOTERM))
df = df[df$GO%in%ref,]
df = merge(df, goterms, by='GO', all.x=TRUE)

vari_df = VariableFeaturePlot(OSN)
vari_df = vari_df$data

res_list = list()
for (i in 1:length(ref)) {
  tmp_GO = df[df$GO==ref[i],]
  tmp_vari_df = vari_df[rownames(vari_df)%in%tmp_GO$SYMBOL,]
  res_list[[i]] = c(tmp_GO$GO_name[1],mean(tmp_vari_df$mean),mean(tmp_vari_df$variance.standardized))
}
results = rbind.data.frame(res_list, make.row.names = F)
results = t(results)
rownames(results) = NULL
colnames(results) = c("GO","mean_expression", "mean_standardized_variance")
results = as.data.frame(results)
results$mean_expression = as.numeric(results$mean_expression)
results$mean_standardized_variance = as.numeric(results$mean_standardized_variance)

ggplot(results, aes(mean_standardized_variance, reorder(GO, mean_standardized_variance)))+geom_bar(stat = 'identity')

# 1D
MOE = readRDS("data/MOE.rds")
Idents(MOE) = MOE$cell_type
MOE_sub = subset(MOE, ident = c("OSN","iOSN"))
VlnPlot(MOE_sub, c("Nrxn1","Nrxn2","Nrxn3"))
VlnPlot(MOE, c("Nrxn1","Nrxn2","Nrxn3"), pt.size = 0)

# 1E
OSN = readRDS("data/OR_OSN_Nrxn.rds")
Idents(OSN) = OSN$OR_identity
OSN_ave = AverageExpression(OSN, return.seurat = T)
mtx = GetAssayData(OSN_ave, slot = "scale.data")[grep("Nrxn",rownames(OSN), value = T),]
mtx[mtx>=2] = 2
mtx[mtx<= (-2)] = (-2)
myColor <- colorRampPalette(c("blue", "white", "Red"))(50)
myBreaks <- c(seq(min(mtx), 0, length.out = 26), seq(0.1, max(mtx), , length.out = 24)) 
pheatmap(mtx, show_colnames = F)

