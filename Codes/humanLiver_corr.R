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


##### preparing the rat dataset #####

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
#old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(new_data_scCLustViz_object)
merged_samples <- your_scRNAseq_data_object
rownames(merged_samples)[(!rownames(merged_samples) %in% mapper$V2)]

merged_samples$res.0.6 <- as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$cluster <- as.character(merged_samples$res.0.6)
table(merged_samples$orig.ident)
cluster_names = merged_samples$cluster 
cluster_names_types = names(table(cluster_names))

DefaultAssay(merged_samples) <- 'RNA'
merged_samples <- FindVariableFeatures(merged_samples)
rat_HVGs <- VariableFeatures(merged_samples)

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
colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
head(cluster_average_exp_df)


#####  Selecting the Hepatocyte populations ##### 
#### final clusters - old samples: resolution 0.6, 17 clusters
### Heps: any clusters other than 13, 9, 5, 7, 14, 3 >> +10, 11 >> final Hep clusters: "0"  "1"  "2"  "4"  "6"  "8"  "12" "15" "16"

#### final clusters - new samples: resolution 0.6, 18 clusters
### Heps:  0, 1, 2, 3, 5, 9, 14, 15 >> 3, 5 not Heps >> final Hep clusters:  0, 1, 2, 9, 14, 15

### QC-refined old samples (unified QC and MT removed)
old_sample_Hep_clusters <- cluster_names_types[!cluster_names_types %in% c('13', '9', '5', '7', '14', '3')] 
clusters_to_check <- paste0('cluster_',old_sample_Hep_clusters)

### MT removed new samples
new_sample_Hep_clusters <-  as.character(c(0, 1, 2, 3, 5, 9, 14, 15))
clusters_to_check <- paste0('cluster_',new_sample_Hep_clusters)

# including only the highly variable features 
#cluster_average_exp_df <- cluster_average_exp_df[rat_HVGs,]
#dim(cluster_average_exp_df)

## scale and center all the genes in the matrix
Hep_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df[,colnames(cluster_average_exp_df) %in% clusters_to_check]) ## scaling over the clusters of interest
Hep_cluster_average_exp$rat_ID=rownames(Hep_cluster_average_exp)
head(Hep_cluster_average_exp)

saveRDS(Hep_cluster_average_exp, 'Results/rat_new_cluster_average_exp_df.rds') #old

#rat_old_cluster_average_exp_df <- readRDS('Results/rat_old_cluster_average_exp_df.rds')
rat_old_cluster_average_exp_df <- readRDS('Results/rat_new_cluster_average_exp_df.rds')

dim(Hep_cluster_average_exp)
head(Hep_cluster_average_exp)




############## generating the human gene expression matrix
#devtools::install_github("BaderLab/HumanLiver")
#library(HumanLiver)
#viewHumanLiver()
HumanLiverSeurat <- readRDS('~/RatLiver/Objects/HumanLiverSeurat.rds')
sCVdL <- readRDS('~/RatLiver/Objects/HumanLiver_supFile.rds')

sample = HumanLiverSeurat
colnames(sample@meta.data)
DefaultAssay(sample) <- 'RNA'

cluster_names = as.character(Idents(sample) )
cluster_names_types = unique(cluster_names)

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  sample[,cluster_names == a_cluster_name]
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, dim)
names(cluster_expression) = paste0('cluster_', names(cluster_expression))
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

write.csv(cluster_average_exp_df, 'Objects//clusterAverageExp_humanLiver_hvg.txt')
#human_cluster_average <- read.csv('~/HumanLiver/pathway/gsva_inputs/clusterAverageExp_humanLiver.txt', header = T)
human_cluster_average <- read.csv('Objects//clusterAverageExp_humanLiver_hvg.txt', header = T)


human_cluster_average.df <- cbind(X=human_cluster_average[,1], 
                                  get_scaled_by_gene(human_cluster_average[,-1]))

human_hep_clusters <- c(5, 14, 6, 15, 1, 3)
human_cluster_average.df <- human_cluster_average.df[,c('X',paste0('cluster_',human_hep_clusters))]
head(human_cluster_average.df)

head(rat_cluster_average.df)
head(human_cluster_average.df)

boxplot(rat_cluster_average.df[,-ncol(rat_cluster_average.df)])
boxplot(human_cluster_average.df[,-1])


## converting the rat genes to human gene names
mapper <- read.table('Data/rat_Lew_01/features.tsv.gz')
head(mapper)
mapper$cleanSymbol <- sapply(strsplit(mapper$V2, '\\.'), '[[', 1)

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
rat_to_human_genes <- rat_to_human_genes[rat_to_human_genes$hsapiens_homolog_associated_gene_name!='',c('symbol', 'hsapiens_homolog_associated_gene_name')]
head(rat_to_human_genes)
human_cluster_average.df.homolog <- merge(human_cluster_average.df, rat_to_human_genes, by.x='X',
                                          by.y='hsapiens_homolog_associated_gene_name', sort=F,all.x=T)
head(human_cluster_average.df.homolog)


#################### 
colnames(human_cluster_average.df.homolog) <- paste0(colnames(human_cluster_average.df.homolog), '.h')
colnames(Hep_cluster_average_exp) <- paste0(colnames(Hep_cluster_average_exp), '.r')

merged_cluster_averages <- merge(human_cluster_average.df.homolog, Hep_cluster_average_exp, 
                                 by.x='symbol.h', by.y='rat_ID.r')
head(merged_cluster_averages)
dim(merged_cluster_averages)

cor.mat <- cor(merged_cluster_averages[,-c(1,2)])
colnames(cor.mat)
paste0('cluster_',c(5, 14, 6, 15, 1, 3),'.h')
pheatmap(cor.mat[paste0('cluster_',c(5, 14, 6, 15, 1, 3),'.h'),7:ncol(cor.mat)], cluster_rows = F)



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





