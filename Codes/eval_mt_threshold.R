#### info from all the input sampels 
folder_name = 'Results/old_samples/MT_info/'
meta_qc_list = lapply(X = list.files(folder_name,pattern = 'meta_*', full.names = T),FUN = readRDS)
names(meta_qc_list) = gsub('meta_qc_df_', '', list.files(folder_name, pattern = 'meta_*'))
meta_qc_list = lapply(meta_qc_list, function(x) {x$id = paste0(x$umi, '_', rownames(x)); x})
lapply(meta_qc_list, dim)
meta_qc_merged = Reduce(rbind,meta_qc_list)

##### merged object info + adding the cluster annotations ##### 
cell_info_df = readRDS('Results/old_samples/MT_info/cell_info_df.rds')
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))


sum(colnames(merged_samples) != cell_info_df$cell_id)
cell_info_df$cluster = merged_samples$cluster
head(cell_info_df)
dim(cell_info_df)
#########################

meta_qc_merged$is_included = meta_qc_merged$id %in% cell_info_df$cell_id
summary(meta_qc_merged$mito_perc[meta_qc_merged$is_included]) #max = 39.96
summary(meta_qc_merged$mito_perc[!meta_qc_merged$is_included]) # max = 66.5

meta_qc_merged_inc = meta_qc_merged[meta_qc_merged$is_included,]
meta_qc_merged_inc2 = merge(cell_info_df, meta_qc_merged_inc, by.x='cell_id', by.y='id')
head(meta_qc_merged_inc2)

df = data.frame(table(paste0('cluster_',meta_qc_merged_inc2$cluster)))
colnames(df) = c('cluster name', 'cell count')
dev.off()
gridExtra::grid.table(df)
#### average mt expression across various clusters
library(dplyr)
group_by(meta_qc_merged_inc2, cluster) %>% summarize(m = mean(mito_perc))

meta_qc_merged_inc_mt30 <- meta_qc_merged_inc2[meta_qc_merged_inc2$mito_perc<30,]
head(meta_qc_merged_inc_mt30)
dim(meta_qc_merged_inc_mt30)
df = data.frame(table(paste0('cluster_',meta_qc_merged_inc_mt30$cluster)))
colnames(df) = c('cluster name', 'cell count')
dev.off()
gridExtra::grid.table(df)


meta_qc_merged_inc_mt20 <- meta_qc_merged_inc2[meta_qc_merged_inc2$mito_perc<20,]
head(meta_qc_merged_inc_mt20)
dim(meta_qc_merged_inc_mt20)
df = data.frame(table(paste0('cluster_',meta_qc_merged_inc_mt20$cluster)))
colnames(df) = c('cluster name', 'cell count')
dev.off()
gridExtra::grid.table(df)


set1_info <- read.csv('figure_panel/set-1-final-info.csv')
gridExtra::grid.table(set1_info[,3])
