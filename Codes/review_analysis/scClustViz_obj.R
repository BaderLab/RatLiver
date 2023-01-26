library(scmap)
library(celldex)
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)


source('~/RatLiver/Codes/Functions.R')
Initialize()


##########################################################################
############## generate scClustViz object for the nuc-seq data ##############
##########################################################################

merged_samples = readRDS('~/rat_sham_sn_data/sham_sn_merged_data.rds')
your_cluster_results = as.character(merged_samples$seurat_clusters)
### calculating the differentially expressed marker genes

## creating the meta.data dataframe ### process was killed
head(merged_samples@meta.data)
your_cluster_results <- data.frame(merged_samples@meta.data[,c(9:ncol(merged_samples@meta.data))]) # 7
rownames(your_cluster_results)  = rownames(merged_samples@meta.data)
head(your_cluster_results)

### calculating the differentially expressed marker genes
sCVdata_list <- CalcAllSCV(
  inD=your_scRNAseq_data_object,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.1, #gene filter - minimum detection rate
  testAll=F, #stop testing clusterings when no DE between clusters
  FDRthresh=0.05,
  calcSil=T, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T
)

#### The altered version used for the Rejection sample
sCVdata_list <- CalcAllSCV(
  inD=merged_samples,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  #exponent=exp(1), #log base of normalized data
  #pseudocount=1,
  #DRthresh=0.5, #gene filter - minimum detection rate
  testAll=T, #stop testing clusterings when no DE between clusters
  #FDRthresh=0.005,
  #calcSil=F, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn= T #
)

#saveRDS(sCVdata_list, '~/LiverTransplant/Objects/sCVdata_list_Integrated_all.rds') ### find the results on run1
#transplant_Integrated_scCLustViz_object <- "~/LiverTransplant/Objects/transplant_Integrated_scCLustViz_object.RData"

transplant_data_scCLustViz_object <- paste0("~/LiverTransplant/Objects/transplant_",sample_name,"_scCLustViz_object.RData")

save(transplant_data,sCVdata_list,
     file=transplant_data_scCLustViz_object) ## new data scClustViz object


load(transplant_data_scCLustViz_object)


runShiny(
  ## write the path to the file bellow:
  filePath= transplant_data_scCLustViz_object,
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  annotationDB="org.Rn.eg.db",
  # This is an optional argument, but will add annotations.
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)


