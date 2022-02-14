source('Codes/Functions.R')
Initialize()



mouse_data_orig <- CreateSeuratObject(counts=Read10X('AcuteLiverFailureDataMouseLiver', gene.column = 2),
                               min.cells=0,
                               min.features=1, 
                               project = "snRNAseq")
sample_order_num = unlist(lapply(str_split(colnames(mouse_data_orig), pattern = '-'), '[[', 2))

mapper = read.table('AcuteLiverFailureDataMouseLiver/genes.tsv')
mouse_annotation = read.table('AcuteLiverFailureDataMouseLiver/cell_annotation.txt')
sample_order = read.csv('AcuteLiverFailureDataMouseLiver/aggregation_sample_order.csv')

mouse_data_orig$label = mouse_annotation$x
sample_2Include = sample_order$order[sample_order$library_id %in% c('SPF_1', 'SPF_2', 'SPF_3')]
mouse_data = mouse_data_orig[,sample_order_num %in% sample_2Include]

head(sample_order)
table(sample_order_num)
head(mouse_annotation)
table(mouse_annotation$x)

# QC based on the paper: cells with <100 detected transcripts and >15% mitochondrial reads were removed
### seems like the filter has already been applied

MIT_PATTERN = '^mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, rownames(mouse_data) )
mouse_data[["mito_perc"]] <- PercentageFeatureSet(mouse_data, features = mito_genes_index)

summary(mouse_data[["mito_perc"]])
summary(mouse_data$nCount_RNA)
summary(mouse_data$nFeature_RNA)


## Normalization
### SCTransform
mouse_data <- SCTransform(mouse_data,conserve.memory=F,
                          verbose=T,return.only.var.genes=F,
                          variable.features.n = nrow(mouse_data), 
                    ## according to the paper scaling is not recommended prior to PCA:
                    ## https://www.biorxiv.org/content/10.1101/576827v2
                    do.scale = FALSE, ### default value 
                    do.center = TRUE) ### default value 

### including only the highly variable genes
mouse_data <- FindVariableFeatures(mouse_data)
mouse_HVGs <- VariableFeatures(mouse_data)
mouse_data = mouse_data[mouse_HVGs,]

######### calculating the cluster average expression values

cluster_names_types = names(table(mouse_data$label))
### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  mouse_data[, mouse_data$label== a_cluster_name]
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, dim)

cluster_average_exp <- lapply(cluster_expression, function(x){
  ## calculate the average expression of each gene in each cluster
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
## Concatenate all the clusters together to make a matrix
mouse_cluster_average_df = do.call(cbind,cluster_average_exp)
#colnames(cluster_average_exp_df) = paste0('cluster_',names(cluster_average_exp))
colnames(mouse_cluster_average_df) = names(cluster_average_exp)
head(mouse_cluster_average_df)
mouse_cluster_average_df$mouse_genes = rownames(mouse_cluster_average_df)
mouse_cluster_average_df <- readRDS('Results/mouse_cluster_average_exp_HVGs.rds')

##### comparison with rat
rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
rat_to_mouse_genes = rat_to_mouse_genes[rat_to_mouse_genes$mmusculus_homolog_orthology_type=='ortholog_one2one',]
rat_to_mouse_genes <- rat_to_mouse_genes[,c('mmusculus_homolog_associated_gene_name', 'symbol')]
dim(rat_to_mouse_genes)
head(rat_to_mouse_genes)

#### adding the rat ortholog gene symbols
mouse_cluster_average.df.homolog <- merge(mouse_cluster_average_df, rat_to_mouse_genes, by.x='mouse_genes',
                                          by.y='mmusculus_homolog_associated_gene_name', sort=F)
head(mouse_cluster_average.df.homolog)

########################################
rat_new_cluster_average_exp_df <- readRDS('Results/rat_new_cluster_average_exp_all.rds')
rat_old_cluster_average_exp_df <- readRDS('Results/rat_old_cluster_average_exp_all.rds')


##### updating the cell-type annotations of set2
colnames_set2 = c("pDC (17)", "Naive T cell (10)", "Erythroid (5)", "Hep (0)" , "Hep (1)", "Hep (14)",
                  "Hep (15)", "Hep (2)", "Hep (3)", "Hep (9)", "Inflammatory Mac (11)", "LSEC (4)",
                  "LSEC (6)", "Non-Inflammatory Mac (8)", "Mature B cell (12)", "gd T cell (7)", 
                  "Non-Inflammatory Mac (13)", "Stellate (16)", "rat_ID" )
colnames(rat_new_cluster_average_exp_df) = colnames_set2


set2_HVGs = readRDS('Results/set2_rat_HVGs.rds')
set1_HVGs = readRDS('Results/set1_rat_HVGs.rds')

merged_set1rat_mouse = merge(rat_old_cluster_average_exp_df[set1_HVGs,], mouse_cluster_average.df.homolog, by.y='symbol', by.x='rat_ID')
merged_set2rat_mouse = merge(rat_new_cluster_average_exp_df[set2_HVGs,], mouse_cluster_average.df.homolog, by.y='symbol', by.x='rat_ID')

merged_set1rat_mouse = merged_set1rat_mouse[,!colnames(merged_set1rat_mouse) %in% c('rat_ID', 'mouse_genes')]
merged_set2rat_mouse = merged_set2rat_mouse[,!colnames(merged_set2rat_mouse) %in% c('rat_ID', 'mouse_genes')]
dim(merged_set1rat_mouse)
dim(merged_set2rat_mouse)

num_set1rat_clusters = ncol(rat_old_cluster_average_exp_df) -1
num_set2rat_clusters = ncol(rat_new_cluster_average_exp_df) -1

#### set-1
colnames(merged_set1rat_mouse)[(num_set1rat_clusters+1):ncol(merged_set1rat_mouse)] = str_sub(string = colnames(merged_set1rat_mouse)[(num_set1rat_clusters+1):ncol(merged_set1rat_mouse)],start = 4)
dim(merged_set1rat_mouse)
cor.mat.set1 <- cor(merged_set1rat_mouse,method = 'pearson')
rownames(cor.mat.set1) <- gsub('Inflammatory', 'Inf', rownames(cor.mat.set1))
pheatmap(cor.mat.set1[1:num_set1rat_clusters,(num_set1rat_clusters+1):ncol(cor.mat.set1)],
         color = inferno(20),main='', clustering_method = 'ward.D') 

#### set-2
colnames(merged_set2rat_mouse)[(num_set2rat_clusters+1):ncol(merged_set2rat_mouse)] = str_sub(string = colnames(merged_set2rat_mouse)[(num_set2rat_clusters+1):ncol(merged_set2rat_mouse)],
                                                                                              start = 4)
dim(merged_set2rat_mouse)
cor.mat.set2 <- cor(merged_set2rat_mouse, method = 'pearson')
pheatmap(cor.mat.set2[1:num_set2rat_clusters,(num_set2rat_clusters+1):ncol(cor.mat.set2)],
         color = inferno(20),main='', clustering_method='ward.D') 

