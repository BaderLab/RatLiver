source('Codes/Functions.R')
source('Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)


# If filter beyond q < 0.05 I'd  -> using an effect size threshold - e.g. (max - min)/min across layers.
# usually end up with 1000 genes from filtering the Halpern data. check if CYP3A4 isn't in your filtered Halpern data then you are filtering too much.
# Another option is to only use genes that are both highly variable in your query data and significant in the Halpern data.

Initialize()

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}


#set1_info <- read.csv('figure_panel/set-1-final-info.csv')
set1_info <- read.csv('figure_panel/set-2-final-info-updated.csv')
colnames(set1_info)[1] = 'clusters'
head(set1_info)

##### preparing the rat dataset #####

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
#old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(new_data_scCLustViz_object)
merged_samples <- your_scRNAseq_data_object
rownames(merged_samples)[(!rownames(merged_samples) %in% mapper$V2)]

merged_samples$res.0.6 <- as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$cluster <- as.character(merged_samples$res.0.6)

merged_samples@meta.data$umi = colnames(merged_samples)
meta.data.label = merge(merged_samples@meta.data, set1_info, by.x='cluster', by.y='clusters')
meta.data.label <- meta.data.label[match(colnames(merged_samples),meta.data.label$umi),]
merged_samples@meta.data = meta.data.label

table(merged_samples$orig.ident)
cluster_names = merged_samples$label 
cluster_names_types = names(table(cluster_names))
cluster_names_types = names(table(merged_samples$label))

DefaultAssay(merged_samples) <- 'RNA'
merged_samples <- FindVariableFeatures(merged_samples)
rat_HVGs <- VariableFeatures(merged_samples)
saveRDS(rat_HVGs, 'Results/set1_rat_HVGs.rds')


### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples, 'data')[,cluster_names == a_cluster_name] 
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)


## calculate the average expression of each gene in each cluster
cluster_average_exp <- lapply(cluster_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(cluster_average_exp, dim)

## Concatenate all the clusters together to make a matrix
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
#colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
colnames(cluster_average_exp_df) = names(cluster_average_exp)
head(cluster_average_exp_df)

cluster_average_exp_df_orig = cluster_average_exp_df

## scale and center all the genes in the matrix
cluster_average_exp_df <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
cluster_average_exp_df$rat_ID=rownames(cluster_average_exp_df)
head(cluster_average_exp_df)

#saveRDS(cluster_average_exp_df, 'Results/rat_new_cluster_average_exp_all.rds') #old
#rat_old_cluster_average_exp_df <- readRDS('Results/rat_old_cluster_average_exp_all.rds')
rat_old_cluster_average_exp_df <- readRDS('Results/rat_new_cluster_average_exp_all.rds')
head(rat_old_cluster_average_exp_df)
rat_cluster_average.df = rat_old_cluster_average_exp_df

colnames(rat_cluster_average.df)
colnames_set2 = c("pDC (17)", "Naive T cell (10)", "Erythroid (5)", "Hep (0)" , "Hep (1)", "Hep (14)",
                  "Hep (15)", "Hep (2)", "Hep (3)", "Hep (9)", "Inflammatory Mac (11)", "LSEC (4)",
                  "LSEC (6)", "Non-Inflammatory Mac (8)", "Mature B cell (12)", "gd T cell (7)", 
                  "Non-Inflammatory Mac (13)", "Stellate (16)", "rat_ID" )

colnames(rat_cluster_average.df) = colnames_set2

############## generating the human gene expression matrix
#devtools::install_github("BaderLab/HumanLiver")
#library(HumanLiver)
#viewHumanLiver()
HumanLiverSeurat <- readRDS('~/RatLiver/Objects/HumanLiverSeurat.rds')
sCVdL <- readRDS('~/RatLiver/Objects/HumanLiver_supFile.rds')

sample = HumanLiverSeurat
colnames(sample@meta.data)
DefaultAssay(sample) <- 'RNA'

### adding the label information to the seurat object
human_label = attr(sCVdL$res.0.8@Clusters, 'ClusterNames')
human_label.df = data.frame(cluster=names(human_label), label=human_label)
head(human_label.df)
human_label.df$label2 = paste0(human_label.df$label, ' (', human_label.df$cluster, ')')

sample@meta.data$umi = colnames(sample)
meta.data.label = merge(sample@meta.data, human_label.df, by.x='res.0.8', by.y='cluster')
meta.data.label <- meta.data.label[match(colnames(sample),meta.data.label$umi),]
sample@meta.data = meta.data.label


cluster_names = as.character(sample$label2)
cluster_names_types = unique(cluster_names)


### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  sample[,cluster_names == a_cluster_name]
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, dim)
names(cluster_expression) = paste0(names(cluster_expression)) #'cluster_', 
cluster_average_exp <- lapply(cluster_expression, function(x){
  ## calculate the average expression of each gene in each cluster
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})

lapply(cluster_average_exp, dim)
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
colnames(cluster_average_exp_df) = names(cluster_average_exp)
head(cluster_average_exp_df)

#### finding the hilghly variable genes
sample <- FindVariableFeatures(sample)
human_HVG <- VariableFeatures(sample)
cluster_average_exp_df <- cluster_average_exp_df[human_HVG,]

#write.csv(cluster_average_exp_df, 'Objects//clusterAverageExp_humanLiver_hvg.txt')
write.csv(cluster_average_exp_df, 'Objects//clusterAverageExp_humanLiver_hvg_allClusters.txt')
#human_cluster_average <- read.csv('~/HumanLiver/pathway/gsva_inputs/clusterAverageExp_humanLiver.txt', header = T)
human_cluster_average <- read.csv('Objects//clusterAverageExp_humanLiver_hvg_allClusters.txt', header = T)



human_cluster_average.df <- cbind(X=human_cluster_average[,1], 
                                  get_scaled_by_gene(human_cluster_average[,-1]))

human_colNames = c("LSECs (12)", "Cholangiocytes (17)", "Macrophages (10)", "ab.T.cells (2)", "Macrophages (4)",
                   "NK.cells (8)", "gd.T.cells (9)", "Hepatocytes (14)","LSECs (13)","gd.T.cells (18)", "LSECs (11)",
                   "Hepatocytes (15)" ,"Mature.B.cells (16)","Hepatic.Stellate.Cells (20)", "Plasma.cells (7)",
                   "Erythroid.cells (19)", "Hepatocytes (3)","Hepatocytes (5)", "Hepatocytes (1)", "Hepatocytes (6)")
colnames(human_cluster_average.df) = c('genes', human_colNames)
head(rat_cluster_average.df)
head(human_cluster_average.df)


boxplot(rat_cluster_average.df[,-ncol(rat_cluster_average.df)])
boxplot(human_cluster_average.df[,-1])

## converting the rat genes to human gene names
#mapper <- read.table('Data/rat_Lew_01/features.tsv.gz')
#head(mapper)
#mapper$cleanSymbol <- sapply(strsplit(mapper$V2, '\\.'), '[[', 1)

# rat_to_human_genes = .getMapped_rat2model_df(ensembl = useDataset('rnorvegicus_gene_ensembl',
#                                                                   mart=useMart("ensembl")),
#                                              candidateGenes = mapper$cleanSymbol,
#                                              model_animal_name = 'hsapiens')
# saveRDS(rat_to_human_genes, '~/RatLiver/Results/rat_to_human_genes_all.rds')

rat_to_human_genes <- readRDS('~/RatLiver/Results/rat_to_human_genes_all.rds')

rat_cluster_average.df$rat_genes = rownames(rat_cluster_average.df)

sum(!rownames(rat_old_cluster_average_exp_df) %in% make.unique(rat_to_human_genes$symbol))/nrow(rat_old_cluster_average_exp_df)

### arounds half of the genes in the human liver map aren't mapped based on rat_to_human_genes
#### around 24% aren't mapped based on the second analysis rat_to_human_genes_all
sum(!human_cluster_average.df$X %in%rat_to_human_genes$hsapiens_homolog_associated_gene_name)/length(human_cluster_average.df$X)

head(rat_to_human_genes)
rat_to_human_genes <- rat_to_human_genes[rat_to_human_genes$hsapiens_homolog_associated_gene_name!='',c('symbol', 
                                                                                                        'hsapiens_homolog_associated_gene_name',
                                                                                                        'hsapiens_homolog_orthology_type')]
head(rat_to_human_genes)
human_cluster_average.df.homolog <- merge(human_cluster_average.df, rat_to_human_genes, by.x='genes',
                                          by.y='hsapiens_homolog_associated_gene_name', sort=F,all.x=T)
head(human_cluster_average.df.homolog)
dim(human_cluster_average.df.homolog)

### only selecting the genes which have one-to-one orthologs
human_cluster_average.df.homolog <- human_cluster_average.df.homolog[human_cluster_average.df.homolog$hsapiens_homolog_orthology_type == 'ortholog_one2one',]
dim(human_cluster_average.df.homolog)
table(human_cluster_average.df.homolog$hsapiens_homolog_orthology_type)

human_cluster_average.df.homolog = human_cluster_average.df.homolog[,-ncol(human_cluster_average.df.homolog)]
dim(human_cluster_average.df.homolog)

#################### 
colnames(human_cluster_average.df.homolog) <- paste0(colnames(human_cluster_average.df.homolog), '.h')
colnames(rat_cluster_average.df) <- paste0(colnames(rat_cluster_average.df), '.r')


#set1_HVGs = readRDS('Results/set1_rat_HVGs.rds')
set2_HVGs = readRDS('Results/set2_rat_HVGs.rds')

rat_cluster_average.df =  rat_cluster_average.df[rat_cluster_average.df$rat_genes.r %in% set2_HVGs,]


merged_cluster_averages <- merge(human_cluster_average.df.homolog, rat_cluster_average.df, 
                                 by.x='symbol.h', by.y='rat_ID.r')
head(merged_cluster_averages)
dim(merged_cluster_averages)

cor.mat <- cor(merged_cluster_averages[,-c(1,2, ncol(merged_cluster_averages))])
#cor.mat = cor.mat[rownames(cor.mat) != 'Cholangiocytes (17).h',colnames(cor.mat) != 'Cholangiocytes (17).h']
colnames(cor.mat)
colnames(cor.mat) = gsub('\\.r', '', gsub(pattern = '\\.h','',colnames(cor.mat)))
colnames(cor.mat) <- gsub('Inflammatory', 'Inf', colnames(cor.mat))
colnames(cor.mat) <- gsub('and', "&\n", colnames(cor.mat))
rownames(cor.mat) = colnames(cor.mat)
#ward.D2 , ward.D2

pheatmap(cor.mat[1:20,21:ncol(cor.mat)],color = inferno(20), clustering_method='ward.D') # the first 20 clusters are for human map



head(merged_cluster_averages)
informative_rat_genes = merged_cluster_averages$rat_genes.r
saveRDS(informative_rat_genes, 'Results/human_set2Rat_corGenes.rds')





######################################################## 
############ rat set-1 and set-2 comparison ############


##### selecting the number of highly variable features to use #####
nFeatures = 500

old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
set1_data <- your_scRNAseq_data_object
DefaultAssay(set1_data) <- 'RNA'
set1_data <- FindVariableFeatures(set1_data, nfeatures=nFeatures)
set1_HVGs <- VariableFeatures(set1_data)
rm(your_scRNAseq_data_object)

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
load(new_data_scCLustViz_object)
set2_data <- your_scRNAseq_data_object
DefaultAssay(set2_data) <- 'RNA'
set2_data <- FindVariableFeatures(set2_data, nfeatures=nFeatures)
set2_HVGs <- VariableFeatures(set2_data)

set1_HVGs == set2_HVGs
#######################################


colnames_set2 = c("pDC (17)", "\u03B1\u03B2 T cell (10)", "Erythroid (5)", "Hep (0)" , "Hep (1)", "Hep (14)",
                  "Hep (15)", "Hep (2)", "Hep (3)", "Hep (9)", "Inflammatory Mac (11)", "LSEC (4)",
                  "LSEC (6)", "Non-Inflammatory Mac (8)", "Mature B cell (12)", "\u03B3\u03B4 T cell (7)", 
                  "Non-Inflammatory Mac (13)", "Stellate (16)", "rat_ID" )

rat_new_cluster_average_exp_df = readRDS('Results/rat_new_cluster_average_exp_all.rds')
colnames(rat_new_cluster_average_exp_df) = colnames_set2
rat_old_cluster_average_exp_df <- readRDS('Results/rat_old_cluster_average_exp_all.rds')

head(rat_new_cluster_average_exp_df)
head(rat_old_cluster_average_exp_df)

#rat_to_human_genes <- readRDS('~/RatLiver/Results/rat_to_human_genes_all.rds')

##### set-2 - filtering genes based on HVGs and one to one orthologs 
rat_new_cluster_average_exp_df <- rat_new_cluster_average_exp_df[rat_new_cluster_average_exp_df$rat_ID %in%  set2_HVGs,]
##### set-1- filtering genes based on HVGs and one to one orthologs
rat_old_cluster_average_exp_df <- rat_old_cluster_average_exp_df[rat_old_cluster_average_exp_df$rat_ID %in% set1_HVGs,]

dim(rat_new_cluster_average_exp_df)
dim(rat_old_cluster_average_exp_df)

colnames(rat_new_cluster_average_exp_df) = paste0(colnames(rat_new_cluster_average_exp_df), '-s2')
colnames(rat_old_cluster_average_exp_df) = paste0(colnames(rat_old_cluster_average_exp_df), '-s1')

rat_merged = merge(rat_old_cluster_average_exp_df, rat_new_cluster_average_exp_df, by.x='rat_ID-s1', by.y='rat_ID-s2')

dim(rat_merged) ## additional layers to ne to be added to the filter -> 10126 one-2-one orthos
rat_cor_mat = cor(rat_merged[,-1],method = 'pearson')
colnames(rat_cor_mat)

colnames(rat_cor_mat) = gsub('\\-s2', '', gsub('\\-s1', '', colnames(rat_cor_mat)))
colnames(rat_cor_mat) <- gsub('Inflammatory', 'Inf', colnames(rat_cor_mat))
rownames(rat_cor_mat) = colnames(rat_cor_mat)
rat_cor_mat = rat_cor_mat[1:17,18:ncol(rat_cor_mat)]
pheatmap(rat_cor_mat,color = inferno(20),  clustering_method='ward.D')


phm =pheatmap(rat_cor_mat,color = inferno(20),  clustering_method='ward.D') # the first 20 clusters are for human map
phmr <- phm$tree_row$order
#### re-ordering two branches in the heatmap rows based on phmr 
reordered_row = c(8,  1,  4,  2,  3,  7,  6,  9, 5, 17, 16, 11, 12, 10, 13, 14, 15)
pheatmap(rat_cor_mat[reordered_row,],color = inferno(20), cluster_rows = F ,clustering_method='ward.D')





#################################
######### Correlation analysis based on the nature commu. supp info data-3 #########

humanCorAnalysisRes <- read.csv('~/RatLiver/humanCorAnalysisRes.csv')
head(humanCorAnalysisRes)

rat_to_human_genes <- readRDS('~/RatLiver/Results/rat_to_human_genes_all.rds')
head(rat_to_human_genes)
sum(!humanCorAnalysisRes$ID %in% rat_to_human_genes$hsapiens_homolog_ensembl_gene)/length(humanCorAnalysisRes$ID)
humanCorAnalysisRes.merged <- merge(humanCorAnalysisRes, rat_to_human_genes, by.x='ID', by.y='hsapiens_homolog_ensembl_gene')
head(humanCorAnalysisRes.merged)
dim(humanCorAnalysisRes.merged)

### only include the one-one prthology groups 
table(humanCorAnalysisRes.merged$hsapiens_homolog_orthology_type)
duplicated(humanCorAnalysisRes.merged$ID)
humanCorAnalysisRes.merged <- humanCorAnalysisRes.merged[humanCorAnalysisRes.merged$hsapiens_homolog_orthology_type=='ortholog_one2one',]

### check if the data is scaled --> it is
boxplot(humanCorAnalysisRes.merged[,paste0('clust',c(1,3,5,6,14,15))])

rat.human.merged <- merge(humanCorAnalysisRes.merged, Hep_cluster_average_exp, by.x='symbol', by.y='rat_ID')

head(Hep_cluster_average_exp)
head(humanCorAnalysisRes.merged)
humanCorAnalysisRes.merged$symbol[!humanCorAnalysisRes.merged$symbol %in% Hep_cluster_average_exp$rat_ID]

human_clusters = c('clust5','clust14','clust6','clust15','clust1','clust3')
### old map
rat_clusters = c('cluster_0','cluster_1','cluster_10','cluster_11','cluster_12', 'cluster_15', 
                 'cluster_16', 'cluster_2', 'cluster_4', 'cluster_6', 'cluster_8')

### new map
rat_clusters =c('cluster_0', 'cluster_1', 'cluster_14', 'cluster_15',
                'cluster_2', 'cluster_3',  'cluster_5', 'cluster_9')


method = 'spearman'
cor_mat <- rcorr(as.matrix(rat.human.merged[,c(human_clusters, rat_clusters)]), type = method)
cor_mat_pVal <- cor_mat$P
cor_mat <- cor_mat$r
cor_mat_pVal.sub <- cor_mat_pVal[human_clusters, rat_clusters]
cor_mat.sub <- cor_mat[human_clusters, rat_clusters]
fdr_mat <- round(cor_mat_pVal.sub, 5)
fdr_mat_char <- ifelse(fdr_mat<0.001, '***', ifelse(fdr_mat<0.01, '**',ifelse(fdr_mat<0.05,'*','') ))

pheatmap(cor_mat.sub,cluster_rows = F, main=paste0('new samples-', method, '-one2one'), #
         fontsize_row = 19,fontsize_col = 19,fontsize_number = 27,
         display_numbers =fdr_mat_char)
head(rat.human.merged)





#############
tmp = merge(rat_old_cluster_average_exp_df, rat_to_human_genes, by.x='rat_ID', by.y='symbol')
one2one_genes = tmp$rat_ID[tmp$hsapiens_homolog_orthology_type == 'ortholog_one2one']
one2one_genes = one2one_genes[!is.na(one2one_genes)]
rat_old_cluster_average_exp_df <- rat_old_cluster_average_exp_df[rat_old_cluster_average_exp_df$rat_ID %in% one2one_genes,] 

