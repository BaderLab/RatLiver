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

merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
ncol1 = ncol(merged_samples@meta.data) 
### perform clustering for a set of resolutions
resolutions = seq(0.4, 2.6, 0.3)
for (res in resolutions){
  merged_samples <- FindClusters(merged_samples, resolution = res, verbose = FALSE)
}
head(merged_samples@meta.data)
your_cluster_results =data.frame(merged_samples@meta.data[,colnames(merged_samples@meta.data) %in% paste0('SCT_snn_res.', resolutions)])
head(your_cluster_results)



### calculating the differentially expressed marker genes

# sCVdata_list <- CalcAllSCV(
#   inD=your_scRNAseq_data_object,
#   clusterDF=your_cluster_results,
#   assayType="SCT", #specify assay slot of data
#   DRforClust="harmony",#reduced dimensions for silhouette calc
#   exponent=exp(1), #log base of normalized data
#   pseudocount=1,
#   DRthresh=0.1, #gene filter - minimum detection rate
#   testAll=F, #stop testing clusterings when no DE between clusters
#   FDRthresh=0.05,
#   calcSil=T, #use cluster::silhouette to calc silhouette widths
#   calcDEvsRest=T,
#   calcDEcombn=T
# )

sCVdata_list <- CalcAllSCV(
  inD=merged_samples,
  clusterDF=your_cluster_results,
  assayType='SCT', #specify assay slot of data
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

#saveRDS(sCVdata_list, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_sCVdata.rds') ### find the results on run1

sham_sn_merged_scCLustViz_object <- paste0("~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC_scCLustViz_object.RData")
#save(merged_samples,sCVdata_list,
#     file=sham_sn_merged_scCLustViz_object) ## new data scClustViz object







load(sham_sn_merged_scCLustViz_object)
runShiny(
  ## write the path to the file bellow:
  filePath= sham_sn_merged_scCLustViz_object,
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  annotationDB="org.Rn.eg.db",
  # This is an optional argument, but will add annotations.
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)


