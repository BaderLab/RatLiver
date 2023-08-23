source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')

library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}

########################################################################
################ Rat - sinle nuc seq data ################
########################################################################
## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')

hep_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_Hep_subclusters_june30.rds')
merged_samples <- hep_data
########################################################################
############### Analysis based on resolution 2.5 ##############

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)


rat_HVGs <- VariableFeatures(merged_samples)


#### defining the mega clusters based on 0.6 resolution
Hep0 = as.character(c(0, 4, 7, 10, 16, 25, 27)) # 27 is the little tail
Hep1 = as.character(c(1, 6, 15, 17, 12, 31)) # 31 is the tail
Hep2 = as.character(c(2, 5, 9, 21)) #
Hep3 = as.character(c(3, 8,13, 23, 32 ))


##########  updated version (V3) - July 19th - FINAL PUBLICATION FIGURE
Hep3 = c("3", "8", "13", "23", "32")
Hep2 = c("21", "9", "2", "5")
Hep0 = c("25", "27", "4", "16", "7", "0", "10")
Hep1 = c("1", "6", "15", "17", "12", "31", "26")


merged_samples$clusters = as.character(merged_samples$SCT_snn_res.2.5)
merged_samples$Hep_clusters = ifelse(merged_samples$clusters %in% as.character(Hep0), 'Hep0', merged_samples$clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep1), 'Hep1', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep2), 'Hep2', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep3), 'Hep3', merged_samples$Hep_clusters) 
table(merged_samples$Hep_clusters)

merged_samples = merged_samples[, merged_samples$Hep_clusters %in% c('Hep0', 'Hep1', 'Hep2', 'Hep3')]

####################################################
##### annotate the hepatocytes based on the mouse zonation layers

merged_samples$cluster <- as.character(merged_samples$Hep_clusters)
#merged_samples$cluster <- as.character(merged_samples$SCT_snn_res.2.5)
cluster_names = merged_samples$cluster 
cluster_names_types = names(table(cluster_names))

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
colnames(cluster_average_exp_df) = paste0(names(cluster_average_exp)) #'cluster_'
head(cluster_average_exp_df)

## scale and center all the genes in the matrix
Hep_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df)
Hep_cluster_average_exp$rat_ID=rownames(Hep_cluster_average_exp)
head(Hep_cluster_average_exp)
dim(Hep_cluster_average_exp)
Hep_cluster_average_exp <- Hep_cluster_average_exp[Hep_cluster_average_exp$rat_ID %in% rat_HVGs,]

rat_cluster_average.df = Hep_cluster_average_exp
head(rat_cluster_average.df)








rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')

########## Importing and cleaning the Halpern dataset ##########
p_value_th = 1e-60 #1e-20#
q_value_th = 1e-25 # 1e-25 worked well 
q_value_th = 0.05
q_value_th = 0.01

q_value_th
### filter the harpen dataset based on q-value
liver_zonation_Halpern_init <- read.csv('~/XSpecies/Data/MouseZonationHalpern/liver_zonation_Halpern.csv')
liver_zonation_Halpern_init$p.values_2 <- ifelse(is.na(liver_zonation_Halpern_init$p.values),1,liver_zonation_Halpern_init$p.values) 
liver_zonation_Halpern_init$q.values_2 <- ifelse(is.na(liver_zonation_Halpern_init$q.values),1,liver_zonation_Halpern_init$q.values) 
mouse_HVGs <- liver_zonation_Halpern_init$Gene.Symbol[ liver_zonation_Halpern_init$q.values_2 < q_value_th]
length(mouse_HVGs)
mouse_HVGs <- unlist(str_split(mouse_HVGs, ';'))

### filtering the genes which are not included in the rat dataset
#liver_zonation_humanMouse <- read.csv('humanLiver_mouseLayers.csv')
included_in_rat_ds <- mouse_HVGs %in% rat_to_mouse_genes$mmusculus_homolog_associated_gene_name
mouse_HVGs = mouse_HVGs[included_in_rat_ds]
length(mouse_HVGs)

#################### 
is_MouseGene_included <- unlist(lapply(str_split(liver_zonation_Halpern_init$Gene.Symbol,';'), function(x) sum(x%in%mouse_HVGs)>0))

liver_zonation_Halpern_filt <- liver_zonation_Halpern_init[is_MouseGene_included,]
liver_zonation_Halpern_filt$Gene.Symbol_cl <- unlist(lapply(str_split(liver_zonation_Halpern_filt$Gene.Symbol,';'), function(x) x[x%in%mouse_HVGs][1]))
head(liver_zonation_Halpern_filt)
dim(liver_zonation_Halpern_filt)

rat_to_mouse_genes <- rat_to_mouse_genes[rat_to_mouse_genes$mmusculus_homolog_associated_gene_name 
                                         %in% liver_zonation_Halpern_filt$Gene.Symbol_cl,]

### adding the gene's meta information
liver_zonation_Halpern_filt <- merge(liver_zonation_Halpern_filt, rat_to_mouse_genes, 
                                     by.x='Gene.Symbol_cl', by.y='mmusculus_homolog_associated_gene_name',all.x=T)
dim(liver_zonation_Halpern_filt)
head(liver_zonation_Halpern_filt) ## the symbol column contains the rat gene names
### scaling the Halpern data
liver_zonation_Halpern_filt_2 <- scale(t(liver_zonation_Halpern_filt[,paste0('Layer.',1:9)]),scale = T,center = T)
liver_zonation_Halpern_filt_2[is.nan(liver_zonation_Halpern_filt_2)] <- 0 
liver_zonation_Halpern_filt_2 <- data.frame(t(liver_zonation_Halpern_filt_2))
liver_zonation_Halpern_filt_2$rat_symbol <- liver_zonation_Halpern_filt$symbol
head(liver_zonation_Halpern_filt_2)
dim(liver_zonation_Halpern_filt_2)


clusters_to_check = paste0('Hep', 0:3)
clusters_to_check = colnames(rat_cluster_average.df)[1:(ncol(rat_cluster_average.df)-1)]
q_value_th

#### making the final merged matrix 
merged_hepExp_mouseLayer=merge(liver_zonation_Halpern_filt_2, Hep_cluster_average_exp, 
                               by.x='rat_symbol',by.y='rat_ID',all.x=F, all.y=F,sort=F)
head(merged_hepExp_mouseLayer)
dim(merged_hepExp_mouseLayer)
sum(duplicated(merged_hepExp_mouseLayer)) ## number of duplicated genes >> probably made by ortholog matching

merged_hepExp_mouseLayer_num <- merged_hepExp_mouseLayer[,colnames(merged_hepExp_mouseLayer) %in% c(paste0('Layer.',1:9),clusters_to_check)]
head(merged_hepExp_mouseLayer_num)
merged_hepExp_mouseLayer_num.m <- reshape2::melt(merged_hepExp_mouseLayer_num)


#### visualizing the results
#pdf('Plots/HalpernCor_newData_25Clusters.pdf')
p0=ggplot2::ggplot(merged_hepExp_mouseLayer_num.m, aes(x=variable, y=value, color=variable))+geom_boxplot()+theme_classic()
print(p0)

#dir_name = '~/rat_sham_sn_data/standardQC_results/Halpern_cor/' 
#dir.create(dir_name)
#saveRDS(merged_hepExp_mouseLayer_num, paste0(dir_name, 'merged_sn_hepExp_mouseLayer_res2.5_damagedHeps_included.rds'))


dim(merged_hepExp_mouseLayer_num)

### calculating correlations
hepExp_mouseLayer_rcorr <- rcorr(as.matrix(merged_hepExp_mouseLayer_num), type="spearman")
#halpern_cor_mat <- cor(merged_hepExp_mouseLayer_num)
halpern_cor_mat <- hepExp_mouseLayer_rcorr$r
halpern_cor_mat_pVal <- hepExp_mouseLayer_rcorr$P

halpern_cor_mat.sub <- halpern_cor_mat[colnames(halpern_cor_mat) %in% clusters_to_check, colnames(halpern_cor_mat) %in% paste0('Layer.',1:9)]
halpern_cor_mat_pVal.sub <- halpern_cor_mat_pVal[colnames(halpern_cor_mat_pVal) %in% clusters_to_check, 
                                                 colnames(halpern_cor_mat_pVal) %in% paste0('Layer.',1:9)]
fdr_mat <- round(halpern_cor_mat_pVal.sub, 5)
fdr_mat_char <- ifelse(fdr_mat<0.001, '***', ifelse(fdr_mat<0.01, '**',ifelse(fdr_mat<0.05,'*','') ))
halpern_cor_mat.sub.t = t(halpern_cor_mat.sub)
colnames(halpern_cor_mat.sub.t) = gsub('er_', ' ', colnames(halpern_cor_mat.sub.t) )
rownames(halpern_cor_mat.sub.t) = gsub('Layer.', 'L', rownames(halpern_cor_mat.sub.t))


halpern_cor_mat.sub.t = halpern_cor_mat.sub.t[,c('Hep1','Hep0','Hep2','Hep3')]

fdr_mat_char.t =t(fdr_mat_char)
fdr_mat_char.t = fdr_mat_char.t[,c('Hep1','Hep0','Hep2','Hep3')]
  
pheatmap::pheatmap(halpern_cor_mat.sub.t, cluster_rows = F, cluster_cols  = F, 
                   display_numbers = fdr_mat_char.t ,
                   #color = plasma(100),
                   fontsize_row = 15,fontsize_col = 15,
                   fontsize_number = 22,
                   fontsize = 12,
                   main='',
                   color=colorRampPalette(c("blue3", "white", "violetred2"))(50)) #inferno(20)


pheatmap(t(df_cor_sub), fontsize =10,fontsize_row=12,fontsize_col=12, main=main, 
         color = colorRampPalette(c("blue3", "white", "red"))(50), )



dev.off()











