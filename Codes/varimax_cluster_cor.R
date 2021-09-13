source('Codes/Functions.R')
Initialize()
###### loading the required libraries #######
library(randomForest)
require(caTools)
library(caret) 
library(tidyr)

########################################
###### loading average expression values
rat_new_cluster_average_exp_df <- readRDS('Results/rat_new_cluster_average_exp_all.rds')
rat_old_cluster_average_exp_df <- readRDS('Results/rat_old_cluster_average_exp_all.rds')

##### filtering based on highly variable genes
set2_HVGs = readRDS('Results/set2_rat_HVGs.rds')
set1_HVGs = readRDS('Results/set1_rat_HVGs.rds')

rat_new_cluster_average_exp_df = rat_new_cluster_average_exp_df[set2_HVGs,]
rat_old_cluster_average_exp_df = rat_old_cluster_average_exp_df[set1_HVGs,]

##### filtering based on top DE genes instead of HVGs which might be too noisy
set1_markers.df = readRDS('Results/old_samples/Cluster_markers_mt40_lib1500_MTremoved.rds')
set1_markers.df <- lapply(set1_markers.df, function(x) x[x$p_val_adj<0.01,])

set2_markers.df = readRDS('Results/new_samples/Cluster_markers_mt_removed.rds')
set2_markers.df <- lapply(set2_markers.df, function(x) x[x$p_val_adj<0.01,])

num_DE_genes = 10
set1_markers = unique(unlist(lapply(set1_markers.df, function(x) rownames(x)[1:num_DE_genes])))
set2_markers = unique(unlist(lapply(set2_markers.df, function(x) rownames(x)[1:num_DE_genes])))

dim(rat_new_cluster_average_exp_df)
rat_old_cluster_average_exp_df = rat_old_cluster_average_exp_df[set1_markers,]
rat_new_cluster_average_exp_df = rat_new_cluster_average_exp_df[set2_markers,]
dim(rat_new_cluster_average_exp_df)

num_PCs = 35
#### loading varimax factors
varimax_res_set1 <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
score_set1 <- data.frame(varimax_res_set1$rotScores)
colnames(score_set1) = paste0('Varimax_', 1:ncol(score_set1))
loading_set1 <- data.frame(varimax_res_set1$rotLoadings[,1:ncol(varimax_res_set1$rotLoadings)])
loading_set1$genes = rownames(loading_set1)

merged_old = merge(rat_old_cluster_average_exp_df,loading_set1[,c(1:num_PCs,ncol(loading_set1))], 
                   by.x='rat_ID', by.y='genes')
dim(merged_old)
merged_old_cor = cor(merged_old[,-1])
cluster_num = ncol(rat_old_cluster_average_exp_df) -1
pheatmap(merged_old_cor[1:cluster_num, (cluster_num+1):ncol(merged_old_cor)],
         cluster_cols = F, 
         main = paste0('set1-Top ',num_DE_genes,' marker cluster-varimax correlation')) #2000HGV



varimax_res_set2 <- readRDS('Results/new_samples/varimax_rotated_object_new.rds')
score_set2 <- data.frame(varimax_res_set2$rotScores)
colnames(score_set2) = paste0('Varimax_', 1:ncol(score_set2))
loading_set2 <- data.frame(varimax_res_set2$rotLoadings[,1:ncol(varimax_res_set2$rotLoadings)])
loading_set2$genes = rownames(loading_set2)


merged_new = merge(rat_new_cluster_average_exp_df,loading_set2[,c(1:num_PCs,ncol(loading_set2))], by.x='rat_ID', by.y='genes')
dim(merged_new)
merged_new_cor = cor(scale(merged_new[,-1]))
cluster_num = ncol(rat_new_cluster_average_exp_df) -1
pheatmap(merged_new_cor[1:cluster_num, (cluster_num+1):ncol(merged_new_cor)],
         cluster_cols = F, 
         main = paste0('set2-Top ',num_DE_genes,' marker cluster-varimax correlation'))



###################################
####### comparing the correlation based approach with the random forest based matching system
result <- readRDS('Objects/cluster_pred_RFs_set1.rds')
result <- readRDS('Objects/cluster_pred_RFs_set2.rds')
RF_models <- result$models
preds <- result$preds
cm_list <- result$cm


matchClust=sapply(RF_models, function(a_model){
  imp.df = data.frame(importance(a_model))        
  imp.df = imp.df[order(imp.df$MeanDecreaseGini, decreasing = T),]
  imp.df$factor = rownames(imp.df)
  colnames(imp.df)[1] = 'cluster'
  imp.df
},simplify = F)


matchClust.df <- do.call(rbind,matchClust)
matchClust.df$clusterName = unlist(lapply(str_split(rownames(matchClust.df),pattern = '\\.'), '[[', 1))
head(matchClust.df)

matchClust.wide <- spread(matchClust.df[,4:6], factor, value = MeanDecreaseGini)
rownames(matchClust.wide) = matchClust.wide[,1]
matchClust.wide = matchClust.wide[,-1]

matchClust.wide = scale(matchClust.wide)

colnames(matchClust.wide) = gsub(pattern = 'var', '', colnames(matchClust.wide))
matchClust.wide = matchClust.wide[,paste0('PC_', 1:50)]

set1_info <- read.csv('figure_panel/set-2-final-info.csv')
set1_info = set1_info[1:18,]
colnames(set1_info)[1] = 'clusters'
matchClust.wide = matchClust.wide[paste0('cluster_',set1_info$clusters),]
rownames(matchClust.wide) = set1_info$label

pheatmap(matchClust.wide, cluster_cols = F)

