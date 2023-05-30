source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)

#### endothelial Halpern zonation analysis - data from: https://www.nature.com/articles/nbt.4231

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



Mesenchymal = c(29, 24) 
Endothelial = c(11, 30) 
clusters_to_include = c(Endothelial, Mesenchymal)

## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')

sham_sn_merged_scCLustViz_object = '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_endothelial_subclusters_sCVdata_res2.5.RData' ### find the results on run1
load(sham_sn_merged_scCLustViz_object)
merged_samples = endo_data
merged_samples$cluster = as.character(merged_samples$SCT_snn_res.0.8)


Resolution = 2.5
resolutions = Resolution
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)

DefaultAssay(merged_samples) <- 'RNA'
merged_samples <- FindVariableFeatures(merged_samples)
rat_HVGs <- VariableFeatures(merged_samples)



##### annotate the hepatocytes based on the mouse zonation layers

cluster_names = as.character(merged_samples$SCT_snn_res.0.8 )
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


clusters_to_check <- paste0('cluster_',clusters_to_include)
clusters_to_check = colnames(cluster_average_exp_df)

## scale and center all the genes in the matrix
endo_cluster_average_exp <- get_scaled_by_gene(cluster_average_exp_df[,colnames(cluster_average_exp_df) %in% clusters_to_check]) ## scaling over the clusters of interest
endo_cluster_average_exp$rat_ID=rownames(endo_cluster_average_exp)
head(endo_cluster_average_exp)

dim(endo_cluster_average_exp)

# including only the highly variable features 
endo_cluster_average_exp <- endo_cluster_average_exp[rat_HVGs,] 

### converting the rat genes to mouse gene names
rat_to_mouse_genes = .getMapped_rat2model_df(ensembl = useDataset('rnorvegicus_gene_ensembl',
                                                                   mart=useMart("ensembl")),
                         candidateGenes = row.names(cluster_average_exp_df),
                         model_animal_name = 'mmusculus')

#rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')

########## Importing and cleaning the Halpern dataset ##########
p_value_th = 1e-60 #1e-20#
p_value_th = 1e-20
q_value_th = 1e-25
q_value_th = 0.05
q_value_th = 0.01
### filter the harpen dataset based on q-value

endo_pcRNA_zones = read.csv('~/rat_sham_sn_data/Endothelial_zonation_Halpern/Endothelial_pcRNAseq_Zonation_table.csv')
endo_scRNA_zones = read.csv('~/rat_sham_sn_data/Endothelial_zonation_Halpern/Endothelial_scRNAseq_Zonation_table.csv')

endo_pcRNA_zones = endo_scRNA_zones
head(endo_pcRNA_zones)

endo_pcRNA_zones$p.values_2 <- ifelse(is.na(endo_pcRNA_zones$Zonation_pvalue),1,endo_pcRNA_zones$Zonation_pvalue) 
endo_pcRNA_zones$q.values_2 <- ifelse(is.na(endo_pcRNA_zones$Zonation_qvalue),1,endo_pcRNA_zones$Zonation_qvalue) 
mouse_HVGs <- endo_pcRNA_zones$Gene_symbol[endo_pcRNA_zones$p.values_2 < q_value_th]
length(mouse_HVGs)

### filtering the genes which are not included in the rat dataset
#liver_zonation_humanMouse <- read.csv('humanLiver_mouseLayers.csv')
included_in_rat_ds <- mouse_HVGs %in% tolower(rat_to_mouse_genes$mmusculus_homolog_associated_gene_name)
mouse_HVGs = mouse_HVGs[included_in_rat_ds]
length(mouse_HVGs)

#################### 
rat_to_mouse_genes$mmusculus_homolog_associated_gene_name = tolower(rat_to_mouse_genes$mmusculus_homolog_associated_gene_name)
### adding the gene's meta information
endo_pcRNA_zones_filt <- merge(endo_pcRNA_zones, rat_to_mouse_genes, 
                                     by.x='Gene_symbol', by.y='mmusculus_homolog_associated_gene_name')

dim(endo_pcRNA_zones)
dim(endo_pcRNA_zones_filt)
head(endo_pcRNA_zones_filt)

endo_pcRNA_zones_filt_2 <- scale(t(endo_pcRNA_zones_filt[,c(paste0('Exp_layer',1:8), paste0('SEM_layer', 1:8))]),scale = T,center = T)
endo_pcRNA_zones_filt_2 <- scale(t(endo_pcRNA_zones_filt[,c(paste0('Exp_zone',1:4), paste0('SEM_zone', 1:4))]),scale = T,center = T)
endo_pcRNA_zones_filt_2[is.nan(endo_pcRNA_zones_filt_2)] <- 0 
endo_pcRNA_zones_filt_2 <- data.frame(t(endo_pcRNA_zones_filt_2))
endo_pcRNA_zones_filt_2$rat_symbol <- endo_pcRNA_zones_filt$symbol
head(endo_pcRNA_zones_filt_2)
dim(endo_pcRNA_zones_filt_2)



#### making the final merged matrix 
merged_endoExp_mouseLayer=merge(endo_pcRNA_zones_filt_2, endo_cluster_average_exp, 
                               by.x='rat_symbol',by.y='rat_ID',all.x=F, all.y=F,sort=F)
head(merged_endoExp_mouseLayer)
dim(merged_endoExp_mouseLayer)
sum(duplicated(merged_endoExp_mouseLayer)) ## number of duplicated genes >> probably made by ortholog matching


merged_endoExp_mouseLayer_num = merged_endoExp_mouseLayer[,-1]
head(merged_endoExp_mouseLayer_num)

library(Hmisc)
### calculating correlations
endoExp_mouseLayer_rcorr <- rcorr(as.matrix(merged_endoExp_mouseLayer_num), type="pearson")
#halpern_cor_mat <- cor(merged_hepExp_mouseLayer_num)
halpern_cor_mat <- endoExp_mouseLayer_rcorr$r
halpern_cor_mat_pVal <- endoExp_mouseLayer_rcorr$P

halpern_cor_mat.sub <- halpern_cor_mat[colnames(halpern_cor_mat) %in% clusters_to_check, !colnames(halpern_cor_mat) %in% clusters_to_check]
halpern_cor_mat_pVal.sub <- halpern_cor_mat_pVal[colnames(halpern_cor_mat_pVal) %in% clusters_to_check, 
                                                 !colnames(halpern_cor_mat) %in% clusters_to_check]
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


markers = c('Lyve1')
i = 1
gene_name = markers[i] 
df_umap <- data.frame(UMAP_1=getEmb(endo_data, 'umap')[,1], 
                      UMAP_2=getEmb(endo_data, 'umap')[,2], 
                      library_size= endo_data$nCount_RNA, 
                      mito_perc=endo_data$mito_perc, 
                      n_expressed=endo_data$nFeature_RNA,
                      cluster=endo_data$SCT_snn_res.0.8, 
                      orig_cluster=endo_data$cluster,
                      cell_status = endo_data$cell_status,
                      nuclear_fraction=endo_data$nuclear_fraction, 
                      Alb=GetAssayData(endo_data)['Alb',], 
                      a_gene = GetAssayData(endo_data)[gene_name,],
                      sample_name = endo_data$sample_name, 
                      strain = endo_data$strain, 
                      umi=colnames(endo_data))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(size=1.5, alpha=0.6)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1.5, alpha=0.6)+theme_classic()

