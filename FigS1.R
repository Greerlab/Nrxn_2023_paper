library(Seurat)
library(ggplot2)
MOE = readRDS("data/MOE.rds")
VlnPlot(MOE, features = c('Nrxn1','Nrxn2','Nrxn3'), group.by = "cell_type")
