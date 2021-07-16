#########################
##### Automation of matching factors and clusters/covariates(strain) ####
##https://www.r-bloggers.com/2018/01/how-to-implement-random-forests-in-r/
source('Codes/Functions.R')
###### loading the required libraries #######
library(randomForest)
require(caTools)
library(caret) 


####### loading the data #######
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
load(new_data_scCLustViz_object)

###########################
##### Training a random forest model to predict the strain labels #####
##########################

varimax_df <- data.frame(readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')$rotScores)
varimax_df <- data.frame(readRDS('~/RatLiver/Results/new_samples/varimax_rotated_object_new.rds')$rotScores)

colnames(varimax_df) <- paste0('varPC_', 1:ncol(varimax_df))
varimax_df$strain <-as.factor(sapply(strsplit(x = rownames(varimax_df), '_'), '[[', 2))
table(varimax_df$strain)

### splitting the train and test data
sample = sample.split(varimax_df$strain, SplitRatio = .75)
train = subset(varimax_df, sample == TRUE) 
test = subset(varimax_df, sample == FALSE)

### training the RF model
strain_pred_RF <- randomForest(strain ~ . , data = train, importance = TRUE)

saveRDS(strain_pred_RF, 'Objects/strain_pred_RF_set2.rds')
#strain_pred_RF <- readRDS('Objects/strain_pred_RF_set1.rds')

pred = predict(strain_pred_RF, newdata=test[,-ncol(test)])
### generating a confusion matrix
cm = table(label=test[,ncol(test)], prediction=pred)
cm
gridExtra::grid.table(cm)
#### evaluating feature importance 
imp.df = data.frame(importance(strain_pred_RF))        
imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
dev.off()
varImpPlot(strain_pred_RF,main = 'Strain Prediction Based on Varimax-PCs')   



###########################
##### Training a random forest model to predict a cluster labels #####
##########################

varimax_df <- data.frame(readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')$rotScores)
varimax_df <- data.frame(readRDS('~/RatLiver/Results/new_samples/varimax_rotated_object_new.rds')$rotScores)

colnames(varimax_df) <- paste0('varPC_', 1:ncol(varimax_df))
cluster.df <- data.frame(cluster=sCVdata_list$res.0.6@Clusters)
cluster.df$cluster <- paste0('cluster_', as.character(cluster.df$cluster))
cluster_list <- levels(as.factor(cluster.df$cluster))

RF_models <- list()
preds <- list()
cm_list = list()

for(aFeature in cluster_list){
  
  varimax_df$feature <- as.factor(ifelse(cluster.df$cluster == aFeature, aFeature, 'other'))
  ### splitting the train and test data
  sample = sample.split(varimax_df$feature, SplitRatio = .75)
  train = subset(varimax_df, sample == TRUE) 
  test = subset(varimax_df, sample == FALSE)
  
  ### training the RF model
  RF_models[[aFeature]] <- randomForest(feature ~ . , data = train, importance = TRUE)
  
  preds[[aFeature]] = predict(RF_models[[aFeature]], newdata=test[,-ncol(test)])
  ### generating a confusion matrix
  cm_list[[aFeature]] = table(test[,ncol(test)], preds[[aFeature]])
  
}

saveRDS(list(models=RF_models, preds=preds, cm=cm_list ), 'Objects/cluster_pred_RFs_set2.rds')
#saveRDS(list(models=RF_models, preds=preds, cm=cm_list ), 'Objects/cluster_pred_RFs_set1.rds')
#result <- readRDS('Objects/cluster_pred_RFs_set1.rds')
result <- readRDS('Objects/cluster_pred_RFs_set2.rds')
RF_models <- result$models
preds <- result$preds
cm_list <- result$cm

#### generating a table from the scores 
accuracy_l <- lapply(cm_list, function(x) {x = confusionMatrix(x);x$overall})
accuracy.df <- data.frame(do.call(cbind,accuracy_l))

byClassScores_l <- lapply(cm_list, function(x) {x = confusionMatrix(x);x$byClass})
byClassScores.df <- data.frame(do.call(cbind,byClassScores_l))

RFclusterEval.df <- rbind(accuracy.df, byClassScores.df)
head(RFclusterEval.df)

#write.csv(RFclusterEval.df, 'Results/old_samples/RFclusterEval_df_set1.csv')
write.csv(RFclusterEval.df, '~/RatLiver/Results/new_samples/RFclusterEval_df_set2.csv')


matchClust=sapply(RF_models, function(a_model){
  imp.df = data.frame(importance(a_model))        
  imp.df = imp.df[order(imp.df$MeanDecreaseGini, decreasing = T),]
  imp.df$factor = rownames(imp.df)
  colnames(imp.df)[1] = 'cluster'
  imp.df[1,]
},simplify = F)
 

matchClust.df <- do.call(rbind,matchClust)
varimax_df$cluster = cluster.df$cluster


pdf('Plots/RF_VarimaxmatchCluster_minGini_set2.pdf', width = 10, height = 10)
for(i in 1:length(RF_models)) { #
  varImpPlot(RF_models[[i]], main=names(RF_models)[i])   
  a_var_factor = matchClust.df$factor[i]
  df = data.frame(Embeddings(your_scRNAseq_data_object, 'umap')[,1:2], 
                  varimax_df, 
                  emb=varimax_df[,a_var_factor])
  df$cluster_num = as.character(sapply(strsplit(df$cluster, '_'), '[[', 2))
  p1=ggplot(df, aes(UMAP_1,UMAP_2, color=emb))+geom_point()+theme_classic()+ggtitle(a_var_factor)+scale_color_viridis(option = 'inferno',direction = -1)
  p2=ggplot(df, aes(UMAP_1,UMAP_2, color=cluster_num))+geom_point()+theme_classic()+ggtitle(a_var_factor)+scale_color_manual(values = colorPalatte)
  p3=ggplot(df, aes(cluster_num,emb, fill=cluster_num))+geom_violin()+theme_classic()+ggtitle(a_var_factor)+scale_fill_manual(values = colorPalatte)
  p4=ggplot(df, aes(varPC_1 ,emb, color=cluster_num))+geom_point()+theme_classic()+ggtitle(a_var_factor)+scale_color_manual(values = colorPalatte)
  gridExtra::grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
  }
dev.off()





