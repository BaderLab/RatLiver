library(scmap)
library(celldex)
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)


source('~/RatLiver/Codes/Functions.R')
Initialize()

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}


merged_samples <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
table(merged_samples$cluster, merged_samples$annot_IM) ## cluster 8 seems to be the macrophage population
mac_cluster_num = 8
##################################################################
######## Subclustering the macrophage cluster

## refined macrophage ids based on subclustering results
macrophage_ids <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_IDs.rds') 

mac_data = merged_samples[, merged_samples$cluster==mac_cluster_num]
#mac_data = merged_samples[, colnames(merged_samples) %in% macrophage_ids]


mac_data_meta = mac_data@meta.data
DefaultAssay(mac_data) <- 'RNA'
mac_data@assays$SCT <- NULL


mac_data = SCTransform(mac_data, vst.flavor = "v2", verbose = TRUE)
mac_data <- RunPCA(mac_data,verbose=T)
plot(100 * mac_data@reductions$pca@stdev^2 / mac_data@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")
mac_data <- RunHarmony(mac_data, group.by.vars = "sample_name", assay.use="RNA")

res = 1.0
mac_data <- RunUMAP(mac_data, reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = res, verbose = FALSE)

dim(mac_data)
markers1 = c('Ptprc','Cd68', 'Cd163', 'Mrc1', 'Clec4f', 'Siglec5')
markers2 = c('Lyz2', 'Clec10a', 'Clec9a') ##'Xcr1' was not captured


markers = c(markers1, markers2)
i = 7
gene_name = markers[i] 
gene_name = 'Lyz2'
df_umap <- data.frame(UMAP_1=getEmb(mac_data, 'umap')[,1], 
                      UMAP_2=getEmb(mac_data, 'umap')[,2], 
                      library_size= mac_data$nCount_RNA, 
                      mito_perc=mac_data$mito_perc, 
                      n_expressed=mac_data$nFeature_RNA,
                      cluster=mac_data$seurat_clusters, 
                      cell_status = mac_data$cell_status,
                      nuclear_fraction=mac_data$nuclear_fraction, 
                      Alb=GetAssayData(mac_data)['Alb',], 
                      a_gene = GetAssayData(mac_data)[gene_name,],
                      sample_name = mac_data$sample_name, 
                      strain = mac_data$strain, 
                      umi=colnames(mac_data))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(size=1.5, alpha=0.6)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle(gene_name)



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=1.7,alpha=0.6)+
  theme_classic()+scale_color_manual(values = c(colorPalatte))+ggtitle(paste0('resolution: ', res))#colorPalatte
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point(size=1.3,alpha=0.7)+theme_classic()

########################################################################################
########### Make the scClustViz object of the macrophage subcluster ##################
########################################################################################

resolutions = seq(0.6, 1.4, 0.2)
for (res in resolutions){
  mac_data <- FindClusters(mac_data, resolution = res, verbose = FALSE)
}


head(mac_data@meta.data)
your_cluster_results =data.frame(mac_data@meta.data[,colnames(mac_data@meta.data) %in% paste0('SCT_snn_res.', resolutions)])
head(your_cluster_results)


sCVdata_list <- CalcAllSCV(
  inD=mac_data,
  clusterDF=your_cluster_results,
  assayType='SCT', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  #exponent=exp(1), #log base of normalized data
  #pseudocount=1,
  #DRthresh=0.5, #gene filter - minimum detection rate
  testAll=T, #stop testing clusterings when no DE between clusters
  #FDRthresh=0.005,
  #calcSil=F, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn= T #
)

sham_sn_merged_scCLustViz_object =  '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_subclusters_sCVdata.RData' ### find the results on run1
#saveRDS(sCVdata_list, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_sCVdata.rds') ### find the results on run1
save(mac_data, sCVdata_list, 
     file=sham_sn_merged_scCLustViz_object) ## new data scClustViz object





######################################################################
#################### sample contribution in each cluster

df <- data.frame(sample_type = mac_data$sample_name, 
                 cluster = as.character(mac_data$SCT_snn_res.1))
rownames(df) = NULL
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")

###### ordering the bars in decreasing order
freq.df = data.frame(table(df$cluster))
freq.df = freq.df[order(freq.df$Freq, decreasing = T),]
cluster_orders = as.character(freq.df$Var1)
counts$cluster= factor(counts$cluster, levels = cluster_orders ) 


ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )

ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Fraction of sample per cell type (%)')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12.5,angle=90,color='black'),
        legend.title = element_blank()) +  
  xlab('')


######################################################################
#############  Calculating the average gene expression ################
######################################################################

merged_samples = mac_data
cluster_names_types = names(table(df_umap$cluster))

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples, 'data')[,merged_samples$SCT_snn_res.1 == a_cluster_name] 
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)

## calculate the average expression of each gene in each cluster
cluster_average_exp <- lapply(cluster_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(cluster_average_exp, dim)
lapply(cluster_average_exp, head)
## Concatenate all the clusters together to make a matrix
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
#colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
colnames(cluster_average_exp_df) = paste0('c_',names(cluster_average_exp))
head(cluster_average_exp_df)

## scale and center all the genes in the matrix
cluster_average_exp_df <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
cluster_average_exp_df$rat_ID=rownames(cluster_average_exp_df)
head(cluster_average_exp_df)

cluster_average_exp_df = cluster_average_exp_df[rownames(cluster_average_exp_df)%in%VariableFeatures(merged_samples),]
head(cluster_average_exp_df)
dim(cluster_average_exp_df)


########################################################################
##################################################################################
###########  Annotation based on similarity with the mouse samples ###############

mouse_cluster_average_df <- readRDS('~/RatLiver/Results/mouse_cluster_average_exp_HVGs.rds')
colnames_to_change = colnames(mouse_cluster_average_df)[1:(ncol(mouse_cluster_average_df)-1)] 
colnames(mouse_cluster_average_df)[1:(ncol(mouse_cluster_average_df)-1)] = substr(colnames_to_change,4,nchar(colnames(mouse_cluster_average_df)))
### converting the IDs
rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
rat_to_mouse_genes = rat_to_mouse_genes[rat_to_mouse_genes$mmusculus_homolog_orthology_type=='ortholog_one2one',]
rat_to_mouse_genes <- rat_to_mouse_genes[,c('mmusculus_homolog_associated_gene_name', 'symbol')]
dim(rat_to_mouse_genes)
head(rat_to_mouse_genes)

#### adding the rat ortholog gene symbols
mouse_cluster_average.df.homolog <- merge(mouse_cluster_average_df, rat_to_mouse_genes, by.x='mouse_genes',
                                          by.y='mmusculus_homolog_associated_gene_name', sort=F)
head(mouse_cluster_average.df.homolog)
head(cluster_average_exp_df) 

merged_rat_mouse = merge(cluster_average_exp_df, mouse_cluster_average.df.homolog, by.y='symbol', by.x='rat_ID')
head(merged_rat_mouse)
dim(merged_rat_mouse)

number_of_clusters = length(cluster_names_types)

cor.mat <- cor(merged_rat_mouse[,!colnames(merged_rat_mouse) %in% c('mouse_genes', 'rat_ID')],method = 'pearson')
rownames(cor.mat) <- gsub('Inflammatory', 'Inf', rownames(cor.mat))
cor.mat = cor.mat[1:number_of_clusters,(number_of_clusters+1):ncol(cor.mat)]
pheatmap(cor.mat,color = inferno(20),main='', clustering_method = 'ward.D2') 



########################################################################
##################################################################
###### importing cluster average means of reference maps 
ref_cluster_average_exp_df <- readRDS('~/RatLiver/Results/rat_new_cluster_average_exp_all.rds')

### refine the annotation names for each of the clusters
colnames_set2 = c("pDC (17)", "Cd3+ (10)", "Erythroid (5)", "Hep (0)" , "Hep (1)", "Hep (14)",
                  "Hep (15)", "Hep (2)", "Hep (3)", "Hep (9)", "Mo/Mac/cDC (11)", "Endo (4)",
                  "Endo (6)", "Marco/Cd5l Mac (8)", "B cell (12)", "gd T cell (7)", 
                  "Marco/Cd5l Mac (13)", "Mesenchymal (16)", "rat_ID" )

colnames(ref_cluster_average_exp_df) = colnames_set2
nFeatures = 2000

new_data_scCLustViz_object <- "~/RatLiver/Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
load(new_data_scCLustViz_object)

ref_data = your_scRNAseq_data_object
DefaultAssay(ref_data) <- 'RNA'
ref_data <- FindVariableFeatures(ref_data, nfeatures=nFeatures)
HVGs <- VariableFeatures(ref_data)

ref_cluster_average_exp_df = ref_cluster_average_exp_df[rownames(ref_cluster_average_exp_df)%in%HVGs,]

####################################
######## merging the Healthy Rat ref data with the new samples
average_exp_merged = merge(cluster_average_exp_df, ref_cluster_average_exp_df, by.x='rat_ID', by.y='rat_ID')

dim(average_exp_merged) ## additional layers to ne to be added to the filter -> 10126 one-2-one orthos
head(average_exp_merged)
rat_cor_mat = cor(average_exp_merged[,-1],method = 'pearson')
colnames(rat_cor_mat)

colnames(rat_cor_mat) <- gsub('Inflammatory', 'Inf', colnames(rat_cor_mat))
rownames(rat_cor_mat) = colnames(rat_cor_mat)
number_of_clusters = length(cluster_names_types)
rat_cor_mat = rat_cor_mat[1:number_of_clusters,(number_of_clusters+1):ncol(rat_cor_mat)] 
pheatmap(rat_cor_mat,color = inferno(20),  clustering_method='ward.D2')


##############################################################################
#####  Finding the DE markers for macrophage subclusters
##############################################################################
Idents(mac_data) = mac_data$SCT_snn_res.1
cluster_names <-  levels(mac_data)
### finding the markers while removing the mito-genes
Cluster_markers <- sapply(1:length(cluster_names), 
                          function(i) FindMarkers(mac_data, 
                                                  ident.1=cluster_names[i],
                                                  logfc.threshold = 0,
                                                  min.pct=0,
                                                  min.cells.feature = 1,
                                                  min.cells.group = 1
                          ), 
                          simplify = FALSE)
names(Cluster_markers) <- cluster_names


P_value_thr = 0.05
Cluster_markers_final <- sapply(1:length(Cluster_markers), function(i) {
  
  ## selecting the cluster of interest's marker dataframe (x)
  x = Cluster_markers[[i]]
  a_cluster_name <- names(Cluster_markers)[i]
  
  ## sort rows of x based on log-fold change
  x = x[order(x$avg_log2FC, decreasing = T),]
  
  ## sort based on adj p-value and sign of logFC  
  # x$ranking_score=-log10(x$p_val_adj+.Machine$double.xmin)*sign(x$avg_log2FC)
  # x = x[order(x$ranking_score, decreasing = T),]
  
  ## add the average expression of each gene as a column
  selected_cells = Idents(mac_data) == a_cluster_name
  data <- GetAssayData(mac_data)[rownames(x),]
  x$avg_exp <- rowSums(data[,selected_cells])/sum(selected_cells)
  
  ## filtering genes with adj-p-value higher than 0.05
  #x = x[x$p_val_adj<P_value_thr,]
  
  return(x)
}, simplify = F)


names(Cluster_markers_final) <- names(Cluster_markers)
lapply(Cluster_markers_final, head)


saveRDS(Cluster_markers_final, paste0('~/rat_sham_sn_data/standardQC_results/sham_sn_macrophage_subcluster_markers_res1.rds'))


#### saving markers #####
dir.create(paste0('~/rat_sham_sn_data/standardQC_results/Mac_subcluster_markers_res1/'))
for(i in 1:length(Cluster_markers_final)){
  #df <- data.frame(genes=rownames(Cluster_markers[[i]]),
  #                 score=Cluster_markers[[i]]$ranking_score)
  df <- data.frame(Cluster_markers_final[[i]])
  print(head(df, 25))
  file_name <- names(Cluster_markers_final)[i]
  write.csv(df, paste0('~/rat_sham_sn_data/standardQC_results/Mac_subcluster_markers_res1/', file_name,'.txt'), 
            col.names = T, row.names = T, quote = F)
}



########################################################################
##############################################################################
#####  Finding the DE genes between the strains within the macrophages
##############################################################################
Idents(mac_data) = mac_data$strain

strain_mac_makers = FindMarkers(mac_data, 
                                ident.1 = 'LEW', 
                                ident.2 = 'DA')

strain_mac_makers$score = -log10(strain_mac_makers$p_val_adj+1e-50)*strain_mac_makers$avg_log2FC
strain_mac_makers_ord = strain_mac_makers[order(strain_mac_makers$avg_log2FC, decreasing = T),]
strain_mac_makers_ord = strain_mac_makers[order(strain_mac_makers$score, decreasing = T),]
strain_mac_makers_ord = strain_mac_makers[order(strain_mac_makers$p_val_adj, decreasing = F),]
head(strain_mac_makers_ord, 30)

##################################################################################
####### Try it with Nebula - if it wasn't working, run a GLMM yourself ###########
library(devtools)
install_github("lhe17/nebula")

library(nebula)
re = nebula(count = input_data2[,1:num_cells_to_include],
            id = designMatrix2$strainLEW,
            pred=designMatrix2[1:num_cells_to_include,][,c(1,2)])

################################
##### Run GSEA on the DE list
strain_mac_makers_ordered = strain_mac_makers[order(strain_mac_makers$score, decreasing = T),]
strain_mac_makers_ordered = data.frame(gene=rownames(strain_mac_makers_ordered), score=strain_mac_makers_ordered$score)
head(strain_mac_makers_ordered,20)
write.table(strain_mac_makers_ordered, '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_strain_DE_genes_MarcoMacs.rnk', 
            quote = F, sep = '\t',row.names = FALSE, col.names = FALSE)

gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'


####### Running GSEA on the markers list ########
rnk_file_path = paste0('~/rat_sham_sn_data/standardQC_results/')
rnk_file <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = T)[1]
working_dir = paste0('~/rat_sham_sn_data/standardQC_results/gsea_results/')
dir.create(working_dir)

print(rnk_file)
analysis_name = gsub('.rnk', '',list.files(rnk_file_path, pattern = '*.rnk')[1])
working_dir = working_dir

GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                     rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                     analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , #set min=15
                     working_dir, " > gsea_output.txt")
system(GSEA_command)


################################################################################################
######## subcluster 0 of the macrophage subclustering result seem to be macrophages ################
mac_data_sub = mac_data[,mac_data$SCT_snn_res.1.2==0]
table(mac_data_sub$sample_name)
macrophage_ids <- colnames(mac_data)[mac_data$SCT_snn_res.1.2==0]
length(macrophage_ids)
macrophage_ids = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_macrophage_IDs.rds')
length(macrophage_ids)



