library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

file1 = '~/Spatial_Rat/MacParland_Catia__Rat_6-1_A1_VIS/'
file2 = '~/Spatial_Rat/MacParland_Catia__Rat_6-1_B1_VIS//'

data = Load10X_Spatial(
  data.dir=file2,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL)

PP1_list= c("Sds", "Hal", "Cyp2f4", "Arg1", "Acly", "Ass1", "Gls2", "Agxt", "Uroc1", "Gldc", "Gls2")
PP2_list= c("Apoc3", "Serpina1", "Apoc1", "Apoe", "Itih4", "Apoa1", "Ttr", "Tf", "Alb", "Pigr", "Orm1", "Rpl3", "Fads1", "Aldh1b1", "Srd5a1", "Hsd17b13")
IZ_list= c("Saa4", "Hint1", "Cyp8b1", "Cyp2e1", "Cox7c", "Fabp1")
CV_list= c("Notum", "Cyp27a1", "Fabp7", "Akr1c1", "Gsta5", "Slc22a1", "Aox3",
           "Sult1e1", "Fmo1", "Oat", "Ahr", "Cyp7a1", "Glul", "Rhbg", "Cyp2e1", "Cyp1a2", "Por")
markers_list = c(CV_list, IZ_list, PP2_list, PP1_list)


pdf('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_A1_VIS_HepMarkers.pdf')
pdf('~/Spatial_Rat/Spatial_MacParland_Catia__Rat_6-1_B1_VIS_HepMarkers.pdf')


for(i in 1:length(markers_list)){
  p=SpatialFeaturePlot(data, features = markers_list[i])+theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
  print(p)
}
dev.off()


plot1 <- VlnPlot(data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(data, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

data <- RunPCA(data, assay = "SCT", verbose = FALSE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, verbose = FALSE)
data <- RunUMAP(data, reduction = "pca", dims = 1:30)
p1 <- DimPlot(data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3)
p1 + p2
SpatialDimPlot(data, cells.highlight = CellsByIdentities(object = data, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)