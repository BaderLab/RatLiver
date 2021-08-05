source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)

############
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
new_data_scCLustViz_object_endothelial <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_EndothelialSub.RData"

#########
load(new_data_scCLustViz_object)
load(old_data_scClustViz_object)

load(new_data_scCLustViz_object_Immune)
load(new_data_scCLustViz_object_endothelial)


merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$cluster = as.character(sCVdata_list$RNA_snn_res.1@Clusters)
merged_samples$sample_name = sapply(str_split(merged_samples$orig.ident, '_'), '[[', 2)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)

merged_samples$EndoSub = ifelse(colnames(merged_samples) %in% 
                                    colnames(your_scRNAseq_data_object), 'Endothelial\ncells', 'other')
#
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      label=merged_samples$EndoSub,
                      Ptprc=GetAssayData(merged_samples)['Ptprc',],
                      Lyve1=GetAssayData(merged_samples)['Lyve1',],
                      Eng=GetAssayData(merged_samples)['Eng',])
                      
#### check if removal of mt-genes would help with the annotation
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      clusters=paste0('', sCVdata_list$res.0.6@Clusters),
                      #subclusters=merged_samples$cluster,
                      sample=merged_samples$sample_name,
                      strain=merged_samples$strain,
                      Ptprc=GetAssayData(merged_samples)['Ptprc',], 
                      umi=colnames(merged_samples))

##### adding the final annotations to the dataframe
set1_info <- read.csv('figure_panel/set-2-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info <- set1_info[1:18,]
set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)


library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point(alpha=0.8, size=1.5)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors)+
  theme(text = element_text(size=15),legend.title = element_blank())#

title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.8, size=3)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors)+
  theme(text = element_text(size=15),legend.title = element_blank())#

nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(8, "Pastel1"))(nb.cols) #Pastel1

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample))+geom_point(alpha=0.8, size=4)+
  theme_classic()+
  #scale_color_manual(values = c("#E69F00", "#009E73"))+
  #scale_color_manual(values = c("pink1", "skyblue1"))+
  #scale_color_brewer(palette = "Set3")+
  scale_color_manual(values=mycolors)+
  theme(text = element_text(size=20),legend.title = element_blank())

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Ptprc))+geom_point(alpha=0.4,size=1.5)+
  scale_color_viridis(direction = -1)+theme_classic()+
  theme(text = element_text(size=20))

##### checking one specific cluster
a_cluster = as.character(c(7, 8, 10, 11, 12, 17)) 
a_cluster = as.character(c(1, 2, 3, 9, 14, 15)) 

df_umap$a_cluster = ifelse(df_umap$clusters %in% a_cluster,  df_umap$clusters, 'other')

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_cluster))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(name='Cluster',values = c('black', 'grey'))+
  theme(text = element_text(size=22))+#,legend.title = element_blank()
  ggtitle(paste0('Cluster-', a_cluster, ' cells over Set-2 UMAP'))


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_cluster))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(name='Clusters',values = c(rep('black', 6), 'grey'))+
  theme(text = element_text(size=22))+#,legend.title = element_blank()
  ggtitle(paste0('Cluster 1, 2, 3, 9, 14, 15 cells\nover Set-2 UMAP'))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(values = c('black', 'grey'))+
  theme(text = element_text(size=22),legend.title = element_blank())+
  ggtitle(paste0('Endothelial subclusters\nSet-2 UMAP'))


################
df <- data.frame(sample_type = merged_samples$sample_name, 
                 cluster = as.character(df_umap$label))
rownames(df) = NULL
head(df)

nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(8, "Pastel1"))(nb.cols)

####### set1 labels
xlabel = c("Hep (0)","Hep (1)","Hep (2)","LSEC (3)","Hep (4)","Non-Inf Mac (5)",
           "Hep (6)","Stellate (7)", "Hep (8)","Inf Mac (9)", "Non-Inf Mac (10)", 
           "LSEC (11)", "Hep (12)", "NK-like & T cells (13)", "Stellate (14)", "Hep (15)", "Hep (16)" )  

set_1_cl_ord = c('Hep (2)'	, 'Hep (4)',  'Hep (8)',  'Hep (16)',  'Hep (1)',  
                 'Hep (6)',  'Hep (0)',  'Hep (12)','Hep (15)',	'LSEC (3)', 
                 'LSEC (11)',	'Stellate (7)','Stellate (14)',	'Non-Inflammatory Mac (5)',		
                 'Non-Inflammatory Mac (10)', 'Inflammatory Mac (9)',	'NK-like and T cells (13)')

xlabel2 = c('Hep (2)'	, 'Hep (4)',  'Hep (8)',  'Hep (16)',  'Hep (1)',  
                 'Hep (6)',  'Hep (0)',  'Hep (12)','Hep (15)',	'LSEC (3)', 
                 'LSEC (11)',	'Stellate (7)','Stellate (14)',	'NonInf Mac(5)',		
                 'NonInf Mac(10)', 'Inf Mac (9)',	'NK-like &\nT cells (13)')

####### set2 labels
set_2_cl_ord = c('Hep (1)', 'Hep (3)', 'Hep (9)', 'Hep (14)', 'Hep (15)', 'Hep (0)', 'Hep (2)', 
  'LSEC (6)', 'LSEC (4)', 'Stellate (16)', 'Mac & DC (8)', 'Non-Inflammatory Mac (13)', 
  'Inflammatory Mac (11)', 'NK cell (7)', 'Cd3 T cell (10)', 'Mature B cell (12)', 
  'B cell (17)', 'Erythroid (5)')

xlabel2 = c('Hep (1)', 'Hep (3)', 'Hep (9)', 'Hep (14)', 'Hep (15)', 'Hep (0)', 'Hep (2)', 
                 'LSEC (6)', 'LSEC (4)', 'Stellate (16)', 'Mac & DC (8)', 'Non-Inf Mac (13)', 
                 'Inf Mac(11)', 'NK cell (7)', 'Cd3 T cell (10)', 'Mature B cell (12)', 
                 'B cell (17)', 'Erythroid (5)')
#### based on sample-type ##### 
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(set1_info$label[1:18]) ) 
counts$cluster= factor(counts$cluster, levels = as.character(set_2_cl_ord) ) 
counts$gcellType = gsub('\\(.*', '', as.character(counts$cluster))
counts$strain=sapply(strsplit(counts$sample_type,'-') ,'[[', 1)
#write.csv(counts, 'figure_panel/set1-counts.csv')


ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  scale_fill_manual(values=mycolors)+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())
  scale_x_discrete(labels = xlabel)+xlab('')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+
  ylab('Fraction of sample per cell type (%)')+xlab('Cluster')+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  #scale_fill_brewer(palette = "Set3")+
  #scale_fill_manual(values=mycolors)+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=11,angle=90,color='black'),
        legend.title = element_blank()) + 
  scale_x_discrete(labels = xlabel2)+xlab('')



########### dodge barplots -->>> why can't I add borders????
#counts2 = data.frame(read.csv('figure_panel/set2-counts.csv'))

ggplot(data=counts, aes(x=gcellType, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        legend.title = element_blank(), axis.text.x = element_text(color='black',size=11, angle=90))+ 
  #scale_x_discrete(labels = c("Hep","Inf Mac", "LSEC", "NK-like\n&T cells", "Non-Inf\nMac","Stellate"))+
  scale_x_discrete(labels = c("B cell ","Cd3\nT cell","Erythroid","Hep","Inf Mac", "LSEC",
                              "Mac\n& DC", "Mature\nB cell","NK cell", "Non-Inf\nMac", "Stellate"))+  
  xlab('')+ylab('Counts')

############## Dot plots ###############
############### Set-1 ########
merged_samples2 = merged_samples

set1_info_ord = set1_info[match(set_2_cl_ord, set1_info$label),]
matrix = as.matrix(set1_info_ord[,4:8])
markers = c()
for(i in 1:nrow(set1_info_ord)) markers = c(markers,matrix[i,])
markers <- markers[markers!='']
markers = unname(markers)
div_num = 53 #64
markers_final = unique(c(markers[1:div_num],'Ptprc', 'Cd68', markers[(div_num+1):length(markers)]))

Idents(merged_samples2) <- factor(df_umap$label, levels = set_2_cl_ord)
DotPlot(merged_samples2, features = markers_final) + RotatedAxis()+ xlab('Markers')+ylab('')+theme(axis.text.x = element_text(size=9))





############## Varimax-related plots ##############
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rot_data <- readRDS('Results/new_samples/varimax_rotated_object_new.rds') ## MT-removed
rot_data <- readRDS('Results/new_samples/immune_varimax_results.rds') ## set2-immune sub-population

rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
embedd_df_rotated <- data.frame(scores)

pc_num = 4
rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                     emb_val=embedd_df_rotated[,pc_num],
                     cluster=as.character(merged_samples$cluster),
                     Library_size=merged_samples$nCount_RNA,
                     num_expressed_genes=merged_samples$nFeature_RNA,
                     strain=merged_samples$strain,
                     sample_name = merged_samples$sample_name)

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=cluster))+geom_point(alpha=0.7,size=4)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(name='Immune\nsubcluster',values = colorPalatte)+theme_classic()+
  theme(text = element_text(size=22))+#, legend.title = element_blank()
  ggtitle(paste0('Set-2 Immune-subclusters over Varimax 1 and ', pc_num))


#### boxplots plots for varimax
ggplot(rot_df, aes(x=cluster, y=emb_val, fill=cluster))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = colorPalatte)+theme_classic()+
  xlab('Immune subcluster')+ylab(paste0('Varimax ', pc_num))+
  theme(text = element_text(size=22), legend.title = element_blank(), legend.position = "none")+
  ggtitle(paste0('Varimax ', pc_num, ' scores over set-2 Immune-subclusters'))


df_umap$Varimax = rot_df$emb_val
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Varimax))+geom_point(alpha=0.5,size=4)+
  scale_color_viridis(direction = -1, option = 'inferno',name=paste0("Varimax ", pc_num))+theme_classic()+
  theme(text = element_text(size=22))+ggtitle(paste0('Varimax ', pc_num, ' scores over set2\nimmune-subclusters UMAP'))
#'plasma', 'inferno', 'magma'



###### strain-specific varimax factors

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.6,size=4.5)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  theme(text = element_text(size=22), legend.title = element_blank())+# "#56B4E9"
  ggtitle(paste0('Distribution of cells based on strain\nover Varimax ', pc_num))
  


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(rot_df, aes(x=strain, y=emb_val, fill=strain))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  xlab('Strain')+ylab(paste0('Varimax ', pc_num))+
  theme(text = element_text(size=22), legend.title = element_blank())+
  stat_summary(fun.data=data_summary)+stat_compare_means(label.x = 0.9, label.y = 10, size=6)

##########################################


###### Kmeans results on set-1 varimax-15 ###### 
kmeans_df <- readRDS('Results/kmeans_var15_1_oldsamples.rds')
ggplot(kmeans_df, aes(x=scores.Varimax_1,y=scores.Varimax_15,color=labels))+geom_point(alpha=0.7,size=4.5)+
  scale_color_manual(values=c("#E69F00","#56B4E9", "#999999"))+theme_classic()+
  theme(text = element_text(size=22), legend.title = element_blank())+# "#56B4E9"
  xlab('Varimax 1')+ylab('Varimax 15')
### check the lables on the total umap
umap_df2 <- data.frame(df_umap,
                       Kmeans_label=as.character(kmeans_df$labels))
ggplot(umap_df2, aes(UMAP_1, UMAP_2,color=Kmeans_label))+geom_point(alpha=0.4,size=3.5)+
  scale_color_manual(values=c("#E69F00","#56B4E9", "#999999"))+theme_classic()+
  theme(text = element_text(size=22), legend.title = element_blank())# "#56B4E9"



#####################################
################# Varimax correlation with technical factors
#####################################
merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

qc_variables <- c( "Library_size", "num_expressed_genes","mito_perc" )
main = 'Varimax-factors correlation with\ntechnical covariates (set-2)'
PCs_to_check <- c(1:25)
embedd_df_rotated_2 <- embedd_df_rotated[,PCs_to_check]
colnames(embedd_df_rotated_2) <- paste0('var-', PCs_to_check)
df.cov <- data.frame(embedd_df_rotated_2,
                     Library_size=merged_samples$nCount_RNA,
                     num_expressed_genes=merged_samples$nFeature_RNA, 
                     mito_perc=merged_samples$mito_perc)
df_cor <- cor(df.cov)
df_cor_sub <- df_cor[colnames(df_cor) %in% qc_variables, !colnames(df_cor)%in% qc_variables]
rownames(df_cor_sub) = c('Library size', '#Expressed\ngenes', 'Mt-gene\nfraction')
pheatmap(t(df_cor_sub), fontsize =10,fontsize_row=12,fontsize_col=12, main=main, 
         color = colorRampPalette(c("navy", "white", "red"))(50), )

pheatmap(t(df_cor_sub), fontsize =10,fontsize_row=12,fontsize_col=12, main=main, 
         color = colorRampPalette(c("blue3", "white", "red"))(50), )




##################################################
############# Random Forest models ################

result <- readRDS('Objects/cluster_pred_RFs_set1.rds')
result <- readRDS('Objects/cluster_pred_RFs_set2.rds')
RF_models <- result$models
preds <- result$preds
cm_list <- result$cm

cluster_num = 13
FeatImpDf <- data.frame(RF_models[[paste0('cluster_', cluster_num)]]$importance)
FeatImpDf$Varimax <-as.character(paste0('var ',1:50))
head(FeatImpDf)

FeatImpDf$Varimax <- factor(FeatImpDf$Varimax, levels=FeatImpDf$Varimax[order(FeatImpDf$MeanDecreaseGini)])
ggplot(FeatImpDf, aes(y=MeanDecreaseGini, x=as.factor(Varimax)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, fill='lightseagreen')+
  coord_flip()+theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 10),  
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15))+xlab('Varimax Factors')+ylab('Mean Decrease Gini Score')+
  ggtitle(paste0('Feature-importance score of a Random-Forest \n model tranied to predict cluster ',cluster_num))


FeatImpDf$Varimax <- factor(FeatImpDf$Varimax, 
                            levels=FeatImpDf$Varimax[order(FeatImpDf$MeanDecreaseGini, decreasing = T)])
FeatImpDf <- FeatImpDf[1:20,]
ggplot(FeatImpDf, aes(y=MeanDecreaseGini, x=as.factor(Varimax)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, fill='lightseagreen')+
  theme_classic()+#coord_flip()
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90),
        axis.text.y = element_text(color = "grey20", size = 10),  
        axis.title.x = element_text(color = "grey20", size = 13),
        axis.title.y = element_text(color = "grey20", size = 13))+
  xlab('Varimax Factors')+ylab('Mean Decrease Gini Score')+
  ggtitle(paste0('Top-20 feature-importance score of a Random-Forest\nmodel tranied to predict cluster ',cluster_num))



strain_pred_RF <- readRDS('Objects/strain_pred_RF_set1.rds')
strain_pred_RF <- readRDS('Objects/strain_pred_RF_set2.rds')


FeatImpDf <- data.frame(strain_pred_RF$importance)
FeatImpDf$Varimax <-as.character(paste0('var ',1:50))
FeatImpDf$Varimax <- factor(FeatImpDf$Varimax, levels=FeatImpDf$Varimax[order(FeatImpDf$MeanDecreaseGini, decreasing = T)])
FeatImpDf <- FeatImpDf[1:20,]
head(FeatImpDf)
ggplot(FeatImpDf, aes(y=MeanDecreaseGini, x=as.factor(Varimax)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, fill='steelblue1')+
  theme_classic()+#coord_flip()
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90),
        axis.text.y = element_text(color = "grey20", size = 10),  
        axis.title.x = element_text(color = "grey20", size = 13),
        axis.title.y = element_text(color = "grey20", size = 13))+xlab('Varimax Factors')+ylab('Mean Decrease Gini Score')+
  ggtitle('Top-20 feature-importance score of a Random-Forest \n model tranied to predict rat strains (set-1)')






############################################################
############### Quality control figures ###############
#############################################
# rat_DA_M09_WK_008_3pr_v3  rat_LEW_M09_WK_009_3pr_v3
# 'rat_DA_M_10WK_003' # 'rat_DA_01_reseq' 'rat_Lew_02'. rat_Lew_01 
INPUT_NAME = 'rat_LEW_M09_WK_009_3pr_v3'
MIT_CUT_OFF = 50
LIB_SIZE_CUT_OFF = 1000
NUM_GENES_DETECTED = 250
# ## import the data
input_from_10x <- paste0("Data/", INPUT_NAME,'/')
seur_raw <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                               min.cells=0,min.features=1, 
                               project = "snRNAseq")

seur_genes_df <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V2, col.name = 'symbol')
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V1, col.name = 'ensembl')

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, seur_raw[['RNA']]@meta.features$symbol )
seur_raw[["mito_perc"]] <- PercentageFeatureSet(seur_raw, features = mito_genes_index)

df = data.frame(library_size= seur_raw$nCount_RNA, 
                mito_perc=seur_raw$mito_perc , 
                n_expressed=seur_raw$nFeature_RNA)


ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point(size=4,alpha=0.7)+
  scale_color_viridis(name='# expressed genes')+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red", size=0.9)+#labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.9)+
  ggtitle('Library size and mitochondrial transcript cutoffs \nSet-2 Lew sample')+
  theme(text = element_text(size=18))


merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2],
                      LibrarySize=merged_samples$nCount_RNA, 
                      numGenes=merged_samples$nFeature_RNA, 
                      MitoExp=merged_samples$mito_perc)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=numGenes))+geom_point(alpha=0.7,size=3.5)+
  scale_color_viridis(name='Library size')+#Mt-gene\nfraction
  theme_classic()+theme(text = element_text(size=17))+
  ggtitle('Library size\nSet-2 map ')#Mitochondrial gene fraction


######################################
################ gProfiler results

gProfiler.df <- read.csv('figure_panel/gProfiler_csv/set2immune-var4-neg300.csv')
gProfiler.df <- gProfiler.df[,c('term_id','negative_log10_of_adjusted_p_value')]
colnames(gProfiler.df) = c('term', 'adj_pVal')
head(gProfiler.df)
dim(gProfiler.df)
gProfiler.df$term <- factor(gProfiler.df$term, levels=gProfiler.df$term[order(gProfiler.df$adj_pVal)])
gProfiler.df <- gProfiler.df[1:30,]

ggplot(gProfiler.df, aes(y=adj_pVal, x=as.factor(term)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, aes(fill=adj_pVal))+# fill='lightseagreen'
  scale_fill_gradient2(name='-log10(p.adj)',low='orange', mid='snow', high='blue')+
  coord_flip()+theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 16),
        axis.text.y = element_text(color = "grey20", size = 16),  
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15))+
  ylab('-log10(adjusted p-value)')+
  xlab(paste0('Enriched TFs based on Varimax-4 top negative genes\n(Set-2 Ptprc+ subpopulation)'))


# set1-var15Pos-chea




#####################################
############## Ucell plots

library(stats)
library(ggpubr)
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))

kmeans_df <- readRDS('Results/kmeans_var15_1_oldsamples.rds')
high_score_cells = rownames(scores)[kmeans_df$labels %in% c('1','2')]


my.matrix <- GetAssayData(your_scRNAseq_data_object)
my.matrix_sub = my.matrix[,colnames(my.matrix)%in%high_score_cells]
gmtFile = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'
geneSets <- getGmt(gmtFile)
names(geneSets)
query_genesets = names(geneSets)[grepl(pattern = 'lymphocyte migration', names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features=geneIds(geneSets[names(geneSets) %in% query_genesets])))
scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
colnames(scores_2)

p1=ggplot(scores_2, aes(y=POSITIVE.REGULATION.OF.LYMPHOCYTE.MIGRATION.GOBP.GO.2000403_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('POSITIVE REGULATION OF LYMPHOCYTE MIGRATION')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"

p2=ggplot(scores_2, aes(y=LYMPHOCYTE.MIGRATION.GOBP.GO.0072676_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('LYMPHOCYTE MIGRATION')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"

p3=ggplot(scores_2, aes(y=REGULATION.OF.LYMPHOCYTE.MIGRATION.GOBP.GO.2000401_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('REGULATION OF LYMPHOCYTE MIGRATION')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
gridExtra::grid.arrange(p1,p2,p3,nrow=1,ncol=3)



query_genesets = names(geneSets)[grepl(pattern = 'response to interferon-beta', names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features=geneIds(geneSets[names(geneSets) %in% query_genesets])))
colnames(scores_2) <- c('CELLULAR.RESPONSE.TO.INTERFERON.BETA', 'RESPONSE.TO.INTERFERON.BETA')
scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])

p1=ggplot(scores_2, aes(y=CELLULAR.RESPONSE.TO.INTERFERON.BETA, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('CELLULAR RESPONSE TO INTERFERON BETA')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
p2=ggplot(scores_2, aes(y=RESPONSE.TO.INTERFERON.BETA, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('RESPONSE TO INTERFERON BETA')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
gridExtra::grid.arrange(p1,p2,nrow=1,ncol=2)



geneSets <- readRDS('~/geneSets_converted2rat.rds') ### CHEA genesets converted to rat genes
a_gene_set_name <- 'SPI1' # 'HNF4A' # GATA1
query_genesets = names(geneSets)[grepl(pattern = a_gene_set_name, names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features= geneSets[names(geneSets) %in% query_genesets]))
colnames(scores_2) <- sapply(strsplit(colnames(scores_2),'\\.'),function(x) paste(x[1],x[5],x[6],sep = '_'))

scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
head(scores_2)

ggplot(scores_2, aes(y=SPI1_HL60_Human_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p1=ggplot(scores_2, aes(y=SPI1_K562_Human_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p2=ggplot(scores_2, aes(y=SPI1_THIOMACROPHAGE_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p3=ggplot(scores_2, aes(y=SPI1_GC_B, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p4=ggplot(scores_2, aes(y=SPI1_ERYTHROLEUKEMIA_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
ggplot(scores_2, aes(y=SPI1_NB4_Human_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
gridExtra::grid.arrange(p1,p2,p4,nrow=1,ncol=3)



a_gene_set_name <- 'STAT4' # 'HNF4A' # GATA1
query_genesets = names(geneSets)[grepl(pattern = a_gene_set_name, names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features= geneSets[names(geneSets) %in% query_genesets]))
colnames(scores_2) <- sapply(strsplit(colnames(scores_2),'\\.'),function(x) paste(x[1],x[5],x[7],sep = '_'))

scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
head(scores_2)

ggplot(scores_2, aes(y=STAT4_TH1_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()




