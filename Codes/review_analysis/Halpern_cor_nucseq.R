source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)

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


## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
Resolution = 2.5
resolutions = Resolution
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)

rat_HVGs <- VariableFeatures(merged_samples)
head(merged_samples)
pheatmap(table(merged_samples$annot_IM_g,merged_samples$SCT_snn_res.2.5))
table(merged_samples$SCT_snn_res.2.5[merged_samples$annot_IM_g == 'Hep'])
nucseq_heps = as.character(c(0:4, 6, 9:11)) ## resolution 0.6

### Hep clusters for resolution 2.5
high_MT_heps = as.character(c(22, 18, 14, 28))
nucseq_heps = c("27", "25", "0", "7", "16", "10", "4", "9", "2", "5", "21","20",
  "23", "3", "8", "13", "32", "26", "12", "17", "7","6", "1", "15", "31")
nucseq_heps = c(nucseq_heps, high_MT_heps)


##### annotate the hepatocytes based on the mouse zonation layers

merged_samples$cluster <- as.character(merged_samples$SCT_snn_res.2.5)



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
colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
head(cluster_average_exp_df)

### QC-refined old samples (unified QC and MT removed)
clusters_to_check <- paste0('cluster_',nucseq_heps)

## scale and center all the genes in the matrix
Hep_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df[,colnames(cluster_average_exp_df) %in% clusters_to_check]) ## scaling over the clusters of interest
Hep_cluster_average_exp$rat_ID=rownames(Hep_cluster_average_exp)
head(Hep_cluster_average_exp)

dim(Hep_cluster_average_exp)

# including only the highly variable features 
# Hep_cluster_average_exp <- Hep_cluster_average_exp[rat_HVGs,] 

## make row names as hugo symbols instead of ensemble
# df <- data.frame(ensembl=rownames(cluster_average_exp_df))
# df <- merge(x=df,y=mapper, by.x='ensembl',by.y='V1', sort=F)
# df <- df[order(match(df$ensembl, rownames(cluster_average_exp_df))),]
# head(df)
# df$ensembl == rownames(cluster_average_exp_df)
# row.names(cluster_average_exp_df) <- make.unique(df$V2)



### converting the rat genes to mouse gene names
# rat_to_mouse_genes = .getMapped_rat2model_df(ensembl = useDataset('rnorvegicus_gene_ensembl',
#                                                                   mart=useMart("ensembl")),
#                         candidateGenes = row.names(cluster_average_exp_df),
#                         model_animal_name = 'mmusculus')

rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')

########## Importing and cleaning the Halpern dataset ##########
p_value_th = 1e-60 #1e-20#
q_value_th = 1e-25
q_value_th = 0.05
q_value_th = 0.01
### filter the harpen dataset based on q-value
liver_zonation_Halpern_init <- read.csv('~/XSpecies/Data/MouseZonationHalpern/liver_zonation_Halpern.csv')
liver_zonation_Halpern_init$p.values_2 <- ifelse(is.na(liver_zonation_Halpern_init$p.values),1,liver_zonation_Halpern_init$p.values) 
liver_zonation_Halpern_init$q.values_2 <- ifelse(is.na(liver_zonation_Halpern_init$q.values),1,liver_zonation_Halpern_init$q.values) 
mouse_HVGs <- liver_zonation_Halpern_init$Gene.Symbol[ liver_zonation_Halpern_init$q.values_2 < q_value_th]
length(mouse_HVGs)
mouse_HVGs <- unlist(str_split(mouse_HVGs, ';'))


check_mouse_orthologs <- c('Cyp2c37', 'Cyp2c50', 'Cyp2c54', 'Cyp3a11', 'Cyp3a16', 'Cyp3a41a',
                           'Cyp3a41b', 'Cyp3a44', 'Cyp3a57', 'Cyp3a59', 'Cyp4a12a', 'Cyp4a12b')
CYP3A4_orth <- c('Cyp3a11', 'Cyp3a16', 'Cyp3a41a', 'Cyp3a41b','Cyp3a44')

check_mouse_orthologs[check_mouse_orthologs %in% mouse_HVGs]
CYP3A4_orth[CYP3A4_orth %in% mouse_HVGs] ## Cyp3a11 is the only ortholog which is present
liver_zonation_Halpern_init$q.values_2[liver_zonation_Halpern_init$Gene.Symbol %in% CYP3A4_orth]

# human CYP2C19 is homologous to three mouse P450 isoforms, namely, Cyp2c37, Cyp2c50, and Cyp2c54; 
# human CYP3A4 is homologous to five mouse P450 isoforms, namely, Cyp3a11, Cyp3a16, Cyp3a41a, Cyp3a41b, and Cyp3a44; 
# human CYP3A43 is homologous to mouse Cyp3a57 and Cyp3a59; 
# human CYP4A22 is homologous to mouse Cyp4a12a and Cyp4a12b 

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



#### making the final merged matrix 
merged_hepExp_mouseLayer=merge(liver_zonation_Halpern_filt_2, Hep_cluster_average_exp, 
                               by.x='rat_symbol',by.y='rat_ID',all.x=F, all.y=F,sort=F)
head(merged_hepExp_mouseLayer)
dim(merged_hepExp_mouseLayer)
sum(duplicated(merged_hepExp_mouseLayer)) ## number of duplicated genes >> probably made by ortholog matching

merged_hepExp_mouseLayer_num <- merged_hepExp_mouseLayer[,colnames(merged_hepExp_mouseLayer) %in% c(paste0('Layer.',1:9),clusters_to_check)]
head(merged_hepExp_mouseLayer_num)
merged_hepExp_mouseLayer_num.m <- melt(merged_hepExp_mouseLayer_num)


#### visualizing the results
#pdf('Plots/HalpernCor_newData_25Clusters.pdf')
p0=ggplot2::ggplot(merged_hepExp_mouseLayer_num.m, aes(x=variable, y=value, color=variable))+geom_boxplot()+theme_classic()
print(p0)

dir_name = '~/rat_sham_sn_data/standardQC_results/Halpern_cor/' 
#dir.create(dir_name)
#saveRDS(merged_hepExp_mouseLayer_num, paste0(dir_name, 'merged_sn_hepExp_mouseLayer_res2.5_damagedHeps_included.rds'))

### calculating correlations
hepExp_mouseLayer_rcorr <- rcorr(as.matrix(merged_hepExp_mouseLayer_num), type="pearson")
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
pheatmap::pheatmap(halpern_cor_mat.sub.t, cluster_rows = F, 
                   display_numbers = t(fdr_mat_char) ,
                   #color = plasma(100),
                   fontsize_row = 15,fontsize_col = 15,fontsize_number = 22,main='',
                   color=colorRampPalette(c("blue3", "white", "violetred2"))(50)) #inferno(20)


pheatmap(t(df_cor_sub), fontsize =10,fontsize_row=12,fontsize_col=12, main=main, 
         color = colorRampPalette(c("blue3", "white", "red"))(50), )



dev.off()












#### checking the UMAP of the integrated samples
plot(100 * merged_sample@reductions$pca@stdev^2 / merged_sample@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

PC_NUMBER = 18
merged_sample <- RunUMAP(merged_sample,dims=1:PC_NUMBER, reduction = "harmony",perplexity=30, reduction.name='umap_h')

umap_emb <- data.frame(Embeddings(merged_samples, 'umap_h'))
umap_emb$sample_type = merged_sample$sample_type
umap_emb$cell_type = merged_sample$cell_type
umap_emb$cluster = as.character(merged_sample$SCT_snn_res.0.3)
colnames(umap_emb)[1:2] = c('UMAP_1', 'UMAP_2')
head(umap_emb)


ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=cell_type),alpha=0.7,size=2)+theme_classic()
ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=cluster),alpha=0.7,size=2)+theme_classic()











###### checking if the important markers are expressed and their pattern in this data
### importing the input markers and removing the 
Total_markers_converted_df <- readRDS('Data/Total_markers_converted_df.rds')
Total_markers_names <- names(Total_markers_converted_df)

### markers which aren't present in the expression matrix
Total_markers_converted_df <- sapply(1:length(Total_markers_converted_df),
                                     function(i){
                                       a_mapped_markers_df <- Total_markers_converted_df[[i]]
                                       a_mapped_markers_df <- a_mapped_markers_df[a_mapped_markers_df$rnorvegicus_homolog_ensembl_gene %in% rownames(merged_sample),]
                                       return(a_mapped_markers_df)
                                     }, simplify = F)

names(Total_markers_converted_df) <- Total_markers_names
plasma_cells_DE[plasma_cells_DE %in% Total_markers_converted_df[['Hepatocytes']]$symbol]
names(Total_markers_converted_df)

a_cluster_merged_markers <- data.frame(merged_markers$cluster_0)

a_cell_type_markers <- Total_markers_converted_df$KCs$rnorvegicus_homolog_ensembl_gene

a_cluster_merged_markers$V2[a_cluster_merged_markers$ensemble_ids %in% a_cell_type_markers]


