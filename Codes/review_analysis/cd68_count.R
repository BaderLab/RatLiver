source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)


old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples@meta.data$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$strain = unlist(lapply(str_split(merged_samples$orig.ident, '_'), '[[', 2))
dim(merged_samples)s
head(merged_samples)
table(merged_samples$cluster)

### adding the cell-type annotation ###
annotation_TLH = read.csv('figure_panel/TLH_annotation.csv')

merged_samples_sub = merged_samples[,merged_samples$cluster %in% c(5, 10, 9)]
df = data.frame(Cd68=GetAssayData(merged_samples_sub)['Cd68',],
                Itgal=GetAssayData(merged_samples_sub)['Itgal',],
                strain=as.character(merged_samples_sub$strain))
df$Cd68_bin = ifelse(df$Cd68 > 0 , 1, 0)
df$Itgal_bin = ifelse(df$Itgal > 0 , 1, 0)
table(df$strain, df$Cd68_bin)
table(df$strain, df$Itgal_bin)
ggplot(df, aes(x=strain, y=Cd68, fill=strain))+geom_boxplot()+theme_classic()
ggplot(df, aes(x=strain, y=Itgal, fill=strain))+geom_boxplot()+theme_classic()


