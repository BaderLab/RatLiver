source('Codes/Functions.R')
Initialize()
library(seriation)

## code for the call-back function
# https://stackoverflow.com/questions/54762405/pheatmap-re-order-leaves-in-dendogram
library(seriation)
cl_cb <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}

#### imputed version #####
# merged_samples_imputed <- readRDS('Results/preproc_rats/merged/merged_samples_imputed.rds')
# merged_samples_norm <- NormalizeData(merged_samples_imputed$estimate)
##########################

merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')
merged_samples_norm <- NormalizeData(merged_samples)
merged_samples_norm <- ScaleData(merged_samples_norm)


#pdf('plots/heatmap_topLoeadingGenes.pdf')
rot_data <- readRDS('Results/preproc_rats/merged/rot_data.rds')
scores <- data.frame(rot_data$rotScores)
Loadings <- rot_data$rotLoadings

num_genes <- nrow(Loadings)
num_cells <- nrow(scores)

colnames(scores) = paste0('PC_', 1:ncol(scores))
scores$UMI <- rownames(scores)

### sort the cells based on their PC embedding > get the indices
PC_value <- 5 #13
num_top_genes <- 30
num_top_cells <- 55

#### generate cell data frame #### 
cell_df <- data.frame(umi=colnames(merged_samples), 
                      PC_1 = scores$PC_1,
                      PC_score = scores[,PC_value],
                      PC_score_rank = rank(scores[,PC_value]))

cell_df$select= cell_df$PC_score_rank %in% c( 1:num_top_cells, (num_cells-num_top_cells):num_cells )
## sanity check
ggplot(cell_df, aes(PC_1, PC_score, color=select))+geom_point(size=4,alpha=0.5)+
  theme_classic()+ylab(paste0('PC_', PC_value))+scale_color_manual(values=c('darkblue', 'red'))+
  theme(legend.title =element_blank(),
        axis.text=element_text(size=16,face="bold"), 
        axis.title=element_text(size=20,face="bold"))
## order the data frame based on the pc-score ranking
cell_df_ord <- cell_df[order(cell_df$PC_score_rank, decreasing = F),]
cell_df_ord <- cell_df_ord[cell_df_ord$select,]

################################################################################

#### generate gene data frame #### 
gene_df <- data.frame(gene=rownames(Loadings), 
                      loading_1 = Loadings[,1],
                      PC_loading=Loadings[,PC_value],
                      PC_loading_rank=rank(Loadings[,PC_value]))

gene_df$select = gene_df$PC_loading_rank %in% c( 1:num_top_genes, (num_genes-num_top_genes):num_genes)

## top-loaded gene in negative side
gene_df[gene_df$PC_loading_rank==nrow(gene_df),]
## top-loaded gene in positive side
gene_df[gene_df$PC_loading_rank==1,]

## order the data frame based on the pc-score ranking
gene_df_ord <- gene_df[order(gene_df$PC_loading_rank, decreasing = F),]
gene_df_ord <- gene_df_ord[gene_df_ord$select,]

#### reorder the original data matrix ####
data <- data.frame(GetAssayData(merged_samples_norm[gene_df_ord$gene, cell_df_ord$umi]))
## imputed data:
#data <- data.frame(merged_samples_norm[gene_df_ord$gene, cell_df_ord$umi])
getHead(data)

data$genes <- rownames(data)
data_G.ord <- data[match(gene_df_ord$gene, data$genes),]
rownames(data) == rownames(data_G.ord)
data_G.ord <- data_G.ord[,-ncol(data_G.ord)]

data.t <- data.frame(t(data_G.ord))
data.t$umi <- rownames(data.t)
data_G.U.ord <- data.t[match(cell_df_ord$umi, data.t$umi),]
rownames(data_G.U.ord) == rownames(data.t)
data_G.U.ord <- data_G.U.ord[,-ncol(data_G.U.ord)]

head(data_G.U.ord)
tail(data_G.U.ord)

### initial check if the heatmap
matrix2Vis <- t(as.matrix(data_G.U.ord))
#colnames(matrix2Vis) <- substr(colnames(matrix2Vis),1,10)
#colnames(matrix2Vis) <- paste0(colnames(matrix2Vis), ' (', cluster_umi_df$cluster,')' )
colnames(matrix2Vis)
rownames(matrix2Vis)
heatmap(matrix2Vis, Rowv = NA, Colv = NA)



## visualizing the expression pattern of selected genes in the PC space

gene_df <- data.frame(gene=rownames(Loadings), 
                      loading_1 = Loadings[,1],
                      PC_loading=Loadings[,PC_value],
                      PC_loading_rank=rank(Loadings[,PC_value]))

genes_to_check <- c('Mrc1', 'Fcgr2a', 'Cd163', 'Cd68', 'C1qa', 'Vcam1', 'Cd106', 'Ccl2', 'Itga1', 'Cd61')
num_genes_to_check = 10

gene_df$select = gene_df$gene %in% genes_to_check
gene_df$select = gene_df$PC_loading_rank %in% c( 1:num_genes_to_check, (num_genes-num_genes_to_check):num_genes)
genes_to_check <- row.names(gene_df[gene_df$select,])

pdf('plots/genesToCheck_PC5.pdf') #_topLoadedGenes
cell_df$strain = merged_samples_norm$strain
p=ggplot(cell_df, aes(PC_1, PC_score, color=strain))+geom_point(size=2,alpha=0.8)+
  theme_classic()
print(p)
for(gene_name in genes_to_check ){
  #'Ly6al'  #'Cyp4a2'
  print(gene_name)
  
  if(! gene_name %in% rownames(merged_samples_norm)) next
  cell_df$gene_express <- GetAssayData(merged_samples_norm)[gene_name,]
  
  p=ggplot(cell_df, aes(PC_1, PC_score, color=gene_express))+geom_point(size=2,alpha=0.8)+
    theme_classic()+scale_color_viridis(direction = -1)+ggtitle(gene_name)
  print(p)
}
dev.off()






###### preparing the annotation data frame  ######

cluster_umi_df <- data.frame(umi=colnames(merged_samples), cluster=merged_samples$merged_clusters)
cluster_umi_df <- cluster_umi_df[cell_df_ord$umi,]
cluster_umi_df <- cluster_umi_df[match(cell_df_ord$umi,cluster_umi_df$umi),]
table(cluster_umi_df$cluster)
rownames(data_G.U.ord) == rownames(cluster_umi_df)

cells_AUC_df <- readRDS('~/cells_AUC_df.rds')
### AUCell scores for the imputed data
#cells_AUC_df <- readRDS("~/XSpecies/Results/preproc_rats/merged/AUCell/PPARG_cells_AUC_df_imputed.rds")
getHead(cells_AUC_df)

gene_set_name <- c('PPARG') # HNF4A
is_a_gene_set <- sapply(str_split(colnames(cells_AUC_df),pattern = '\\.'), function(x) x[1]) %in% gene_set_name
matched_gene_set_names <- colnames(cells_AUC_df)[is_a_gene_set]
cells_AUC_df_sub <- cells_AUC_df[rownames(cells_AUC_df) %in% colnames(matrix2Vis),matched_gene_set_names]
cells_AUC_df_sub_ord <- cells_AUC_df_sub[match(colnames(matrix2Vis), rownames(cells_AUC_df_sub)),] 


annotation_df <- data.frame(strain=ifelse(substr(colnames(matrix2Vis),1,10) %in% 
                                            c('rat_Lew_01', 'rat_Lew_02'), 'LEW rat', 'DA rat' ), 
                            cluster=cluster_umi_df$cluster, 
                            cells_AUC_df_sub_ord)

row.names(annotation_df) == rownames(cells_AUC_df_sub_ord)

annotation_df <- annotation_df[,c('strain', 'cluster', 'PPARG.20176806.ChIP.Seq.MACROPHAGES.Mouse')]
colnames(annotation_df)[3] <- 'PPARG.ChIP.Seq.MACROPHAGES.Mouse'
colnames(annotation_df)[3:8]

### checking the distribution of the gene-set's scores in the two strains ####
pdf(paste0('plots/PC',PC_value,'_aucellScores_boxplot2.pdf'), width = 6, height = 8) #_imputed
ggplot(annotation_df, aes(y=PPARG.19300518.ChIP.PET.3T3.L1.Mouse, x=strain, fill=strain))+
  geom_boxplot()+theme_classic()+ggtitle(paste0('PC-', PC_value))
ggplot(annotation_df, aes(y=PPARG.20176806.ChIP.Seq.THIOMACROPHAGE.Mouse, x=strain, fill=strain))+
  geom_boxplot()+theme_classic()+ggtitle(paste0('PC-', PC_value))
ggplot(annotation_df, aes(y=PPARG.23326641.ChIP.Seq.C3H10T1.2.Mouse, x=strain, fill=strain))+
  geom_boxplot()+theme_classic()+ggtitle(paste0('PC-', PC_value))
ggplot(annotation_df, aes(y=PPARG.20176806.ChIP.Seq.3T3.L1.Mouse, x=strain, fill=strain))+
  geom_boxplot()+theme_classic()+ggtitle(paste0('PC-', PC_value))
ggplot(annotation_df, aes(y=PPARG.20176806.ChIP.Seq.MACROPHAGES.Mouse, x=strain, fill=strain))+
  geom_boxplot()+theme_classic()+ggtitle(paste0('PC-', PC_value))
ggplot(annotation_df, aes(y=PPARG.20887899.ChIP.Seq.3T3.L1.Mouse, x=strain, fill=strain))+
  geom_boxplot()+theme_classic()+ggtitle(paste0('PC-', PC_value))
dev.off()



row_font_size = 11
pdf(paste0('plots/PC',PC_value,'_pheatmaps_cellScore2_selectedGenes.pdf'), width = 20, height = 13) # width = 10 height = 12
pheatmap::pheatmap(matrix2Vis, kmeans_k=NA, cluster_rows=F, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=F, show_rownames = T, show_colnames =F,
                   annotation_col = annotation_df,
                   fontsize_row = row_font_size, fontsize_col = 8, border_color='black')

pheatmap::pheatmap(matrix2Vis, kmeans_k=NA, cluster_rows=T, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=F, show_rownames = T, show_colnames =F,
                   annotation_col = annotation_df,
                   fontsize_row = row_font_size, fontsize_col = 8, border_color='black')

### cluster the cells and genes
pheatmap::pheatmap(matrix2Vis, kmeans_k=NA, cluster_rows=T, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=T, show_rownames = T, show_colnames =F,
                   fontsize_row = row_font_size, fontsize_col = 8,
                   annotation_col = annotation_df,
                   clustering_callback = cl_cb)
dev.off()






#### check the expression of enriched gene sets original data matrix ####

PC_5_hits <- c('PPARG', 'CLOCK', 'HNF4A', 'ESR1', 'TAL1', 'EGR1') 
PC_13_hits <- c('GATA1','SPI1' ,'MECOM', 'MITF', 'RUNX1', 'PPARG') #

dir.create(paste0('plots/PC_geneSets_heatmaps/PC',PC_value,'/'))

PC_hits <- PC_5_hits
for(geneSet_Name in PC_hits){
  
  print(geneSet_Name)
  geneSets_converted2rat <- readRDS('~/geneSets_converted2rat.rds')
  ## subsetting the all genes-set of interest
  is_a_geneset <- sapply(str_split(names(geneSets_converted2rat),pattern = ' '), 
                              function(x) x[1]) == geneSet_Name
  geneSets_converted2rat <- geneSets_converted2rat[is_a_geneset]

  
  
  dir.create(paste0('plots/PC_geneSets_heatmaps/PC',PC_value, '/',geneSet_Name,'/'))
  
  for(gene_set_number in 1:length(geneSets_converted2rat)){
    
    print(gene_set_number)
    ### select a single gene-set form the list of gene-sets matched to a term
    a_gene_set <- geneSets_converted2rat[[gene_set_number]]
    a_gene_set_name <- names(geneSets_converted2rat)[gene_set_number]
    
    ## order the data frame based on the pc-score ranking
    gene_df_ord2 <- gene_df[order(gene_df$PC_loading_rank, decreasing = F),]
    gene_df_ord2$select <- gene_df_ord2$gene %in% a_gene_set
    gene_df_ord2$HVG <- gene_df_ord2$gene %in% VariableFeatures(merged_samples)
    gene_df_ord2 <- gene_df_ord2[gene_df_ord2$select ,]
    dim(gene_df_ord2)
    
    if(nrow(gene_df_ord2)>100) gene_df_ord2 <- gene_df_ord2[gene_df_ord2$HVG,]
    dim(gene_df_ord2)
    
    data <- data.frame(GetAssayData(merged_samples_norm[gene_df_ord2$gene, cell_df_ord$umi]))
    
    ### order on genes
    data$genes <- rownames(data)
    data_G.ord <- data[match(gene_df_ord2$gene, data$genes),]
    data_G.ord <- data_G.ord[,-ncol(data_G.ord)]
    ### order on cells
    data.t <- data.frame(t(data_G.ord))
    data.t$umi <- rownames(data.t)
    data_G.U.ord <- data.t[match(cell_df_ord$umi, data.t$umi),]
    data_G.U.ord <- data_G.U.ord[,-ncol(data_G.U.ord)]
    
    ### convert to numeric matrix 
    matrix2Vis <- t(as.matrix(data_G.U.ord))
    
    ### find the AUCell score
    is_a_gene_set <- sapply(str_split(colnames(cells_AUC_df),pattern = '\\.'), function(x) x[1]) %in% geneSet_Name
    matched_gene_set_names <- colnames(cells_AUC_df)[is_a_gene_set]
    cells_AUC_df_sub <- data.frame(cells_AUC_df[rownames(cells_AUC_df) %in% colnames(matrix2Vis),matched_gene_set_names])
    
    cells_AUC_df_sub_rowNames <- rownames(cells_AUC_df)[rownames(cells_AUC_df) %in% colnames(matrix2Vis)]
    cells_AUC_df_sub_ord <- cells_AUC_df_sub[match(colnames(matrix2Vis), cells_AUC_df_sub_rowNames),] 
    
    
    annotation_df <- data.frame(strain=ifelse(substr(colnames(matrix2Vis),1,10) %in% 
                                                c('rat_Lew_01', 'rat_Lew_02'), 'LEW rat', 'DA rat' ), 
                                cluster=cluster_umi_df$cluster, 
                                cells_AUC_df_sub_ord)
    
    row.names(annotation_df) =colnames(matrix2Vis)
    
    width =  22
    height = round(0.03*length(a_gene_set),2)
    pdf(paste0('plots/PC_geneSets_heatmaps/PC',PC_value, '/',geneSet_Name,'/PC',PC_value,'_',a_gene_set_name,'_geneSets_heatmaps.pdf'), 
        width = width, height = height)
    
    pheatmap::pheatmap(matrix2Vis, kmeans_k=NA, cluster_rows=F, 
                       clustering_distance_cols = "manhattan",
                       cluster_cols=F, show_rownames = T, show_colnames =F,
                       fontsize_row = 8.7, fontsize_col = 4,
                       annotation_col = annotation_df,
                       clustering_callback = cl_cb, main = a_gene_set_name)
    
    pheatmap::pheatmap(matrix2Vis, kmeans_k=NA, cluster_rows=T, 
                       clustering_distance_cols = "manhattan",
                       cluster_cols=F, show_rownames = T, show_colnames =F,
                       fontsize_row = 8.7, fontsize_col = 4,
                       annotation_col = annotation_df,
                       clustering_callback = cl_cb, main = a_gene_set_name)
     
    pheatmap::pheatmap(matrix2Vis, kmeans_k=NA, cluster_rows=T, 
                       clustering_distance_cols = "manhattan",
                       cluster_cols=T, show_rownames = T, show_colnames =F,
                       fontsize_row = 8.7, fontsize_col = 4,
                       annotation_col = annotation_df,
                       clustering_callback = cl_cb, main = a_gene_set_name)
    
    dev.off()
  }
}





# scp -r delaram@192.168.142.161:/home/delaram/XSpecies/plots/PC_geneSets_heatmaps/ .



########### heatmap based on AUCell scores
cells_AUC_df <- readRDS('~/cells_AUC_df.rds')
cells_AUC_df[1:3,1:3]
class(cells_AUC_df)

#### reorder the original data matrix ####
data <- data.frame(t(cells_AUC_df)[,cell_df_ord$umi])
data_G.ord <- data
class(data_G.ord)

data.t <- data.frame(t(data_G.ord))
getHead(data.t)
data.t$umi <- rownames(data.t)
data_G.U.ord <- data.t[match(cell_df_ord$umi, data.t$umi),]
rownames(data_G.U.ord) == rownames(data.t)
getHead(data_G.U.ord)
data_G.U.ord <- data_G.U.ord[,!colnames(data_G.U.ord) %in% c('umi', 'UMI')]

head(data_G.U.ord)
colnames(data_G.U.ord)

matrix2Vis <- as.matrix((t(data_G.U.ord)))
matrix2Vis <- data.matrix(data_G.U.ord)

pathway_names <- sapply(colnames(matrix2Vis), function(x)str_split(x, pattern = '\\.')[1])
pathway_names <- sapply(pathway_names, function(x) x[1])
colnames(matrix2Vis) <- pathway_names

class(matrix2Vis)
getHead(matrix2Vis)
#colnames(matrix2Vis) <- substr(colnames(matrix2Vis),1,10)
#colnames(matrix2Vis) <- paste0(colnames(matrix2Vis), ' (', cluster_umi_df$cluster,')' )
heatmap.2(t(matrix2Vis[,1:50]))
heatmap(matrix2Vis, Rowv = NA, Colv = NA)

annotation_df <- data.frame(strain=ifelse(substr(row.names(matrix2Vis),1,10) %in% 
                                            c('rat_Lew_01', 'rat_Lew_02'), 'LEW rat', 'DA rat' ), 
                            cluster=cluster_umi_df$cluster)

row.names(annotation_df) <- row.names(matrix2Vis)
head(annotation_df)


## code for the call-back function
# https://stackoverflow.com/questions/54762405/pheatmap-re-order-leaves-in-dendogram
library(seriation)
cl_cb <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}

row_font_size = 3
pdf(paste0('plots/PC',PC_value,'_pheatmaps_pathway.pdf'), height = 30, width = 8) #width = 10, height = 12
pheatmap::pheatmap(t(matrix2Vis), kmeans_k=NA, cluster_rows=F, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=F,  show_rownames = T, show_colnames =F,
                   annotation_col = annotation_df,
                   fontsize_row = row_font_size, fontsize_col = 8)

### cluster the cells
pheatmap::pheatmap(t(matrix2Vis), kmeans_k=NA, cluster_rows=F, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=T, show_rownames = T, show_colnames =F,
                   fontsize_row = row_font_size, fontsize_col = 8,
                   annotation_col = annotation_df,
                   clustering_callback = cl_cb)
### cluster the genes
pheatmap::pheatmap(t(matrix2Vis), kmeans_k=NA, cluster_rows=T, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=F, show_rownames = T, show_colnames =F,
                   fontsize_row = row_font_size, fontsize_col = 8,
                   annotation_col = annotation_df,
                   clustering_callback = cl_cb)


### cluster both genes and cells
pheatmap::pheatmap(t(matrix2Vis), kmeans_k=NA, cluster_rows=T, 
                   clustering_distance_cols = "manhattan",
                   cluster_cols=T, show_rownames = T, show_colnames =F,
                   fontsize_row = row_font_size, fontsize_col = 8,
                   annotation_col = annotation_df,
                   clustering_callback = cl_cb)


dev.off()



pdf('plots/rat_sample_count_distributions.pdf')
#### Checking if the pattern is not caused by technical artifacts ####
check_df = data.frame(counts=colSums(merged_samples),
           norm_counts= colSums(merged_samples_norm),
           sample=merged_samples$sample_type,
           strain=merged_samples$strain)

ggplot(check_df, aes(counts, fill=strain))+geom_histogram(color="black",alpha=0.4,bins = 40)+
  theme_classic()+ggtitle('Raw counts histogram')+xlab('Sum counts')
ggplot(check_df, aes(counts, fill=sample))+geom_histogram(color="black",alpha=0.4,bins = 40)+
  theme_classic()+ggtitle('Raw counts histogram')+xlab('Sum counts')

ggplot(check_df, aes(y=counts, x=sample,fill=sample))+geom_boxplot(alpha=0.4)+theme_classic()+ggtitle('Raw counts distribution')+xlab('Sum counts')
ggplot(check_df, aes(y=counts, x=strain,fill=strain))+geom_boxplot(alpha=0.4)+theme_classic()+ggtitle('Raw counts distribution')+xlab('Sum counts')


ggplot(check_df, aes(norm_counts, fill=strain))+geom_histogram(color="black",alpha=0.4,bins = 40)+
  theme_classic()+ggtitle('Normalized counts histogram')+xlab('Sum counts')
ggplot(check_df, aes(norm_counts, fill=sample))+geom_histogram(color="black",alpha=0.4,bins = 40)+
  theme_classic()+ggtitle('Normalized counts histogram')+xlab('Sum counts')


ggplot(check_df, aes(y=norm_counts, x=sample,fill=sample))+geom_boxplot(alpha=0.4)+theme_classic()+ggtitle('Normalized counts distribution')+xlab('Sum counts')
ggplot(check_df, aes(y=norm_counts, x=strain,fill=strain))+geom_boxplot(alpha=0.4)+theme_classic()+ggtitle('Normalized counts distribution')+xlab('Sum counts')

dev.off()













