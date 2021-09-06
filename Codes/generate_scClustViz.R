source('Codes/Functions.R')
Initialize()

### importing data and changing the features to gene-symbols
mapper = read.delim('Data/rat_DA_M09_WK_008_3pr_v3/features.tsv.gz',header=F)
#merged_samples <- readRDS('~/RatLiver/Objects/merged_samples_newSamples.rds')
#merged_samples <- readRDS('Objects/merged_samples_oldSamples.rds')

merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

merged_samples <- readRDS('Results/new_samples/Immune_subclusters.rds') 
merged_samples <- readRDS('Results/new_samples/endothelial_subclusters.rds') 
PC_NUMBER = 15

your_scRNAseq_data_object = seur
merged_samples = your_scRNAseq_data_object
assay_data <- GetAssayData(merged_samples, 'data')
#rownames(assay_data) <- mapper$V2



###### new sample (set-2) strain markers werer flipped - fixing the bug #####
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
### flipping the sample tag
merged_samples$strain = ifelse(merged_samples$strain == 'DA', 'LEW', 'DA')
umi_only = sapply(str_split(colnames(merged_samples), '_'), '[[', 6)
umi_sample_info = ifelse(merged_samples$strain == "DA", 'rat_DA_M09_WK_008_', 'rat_LEW_M09_WK_009_')

corrected_UMIs = paste0(umi_sample_info, umi_only)
colnames(assay_data) = corrected_UMIs

#### quality check
colnames(assay_data) == colnames(merged_samples)
colnames(assay_data) == corrected_UMIs




##### removing the mitochondrial genes to eliminate their effect on clustering #####
MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, row.names(merged_samples) )
merged_samples_original <- merged_samples
merged_samples <- merged_samples[-mito_genes_index,]


merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples))  
merged_samples <- RunHarmony(merged_samples, "sample_name",assay.use="RNA")

######## Run this section only for the old data ########
### running UMAP and tSNE on the Old merged samples 
ElbowPlot(merged_samples, 'pca')
PC_NUMBER = 10
merged_samples <- RunUMAP(merged_samples, dims=1:PC_NUMBER, reduction="harmony")
merged_samples <- RunTSNE(merged_samples, dims=1:PC_NUMBER, reduction="harmony")

### evaluate if removing mito-genes will change the umap clustering distribution 
merged_samples <- RunHarmony(merged_samples, "sample_name",assay.use="RNA")
ElbowPlot(merged_samples, 'harmony')
merged_samples <- RunUMAP(merged_samples, dims=1:PC_NUMBER, reduction="harmony")
plot(Embeddings(merged_samples, 'umap')[,1], Embeddings(merged_samples, 'umap')[,2])


merged_samples <- FindNeighbors(merged_samples,
                                reduction="harmony",
                                dims=1:PC_NUMBER,
                                verbose=T)

sCVdata_list <- get_cluster_object(seur = merged_samples, 
                                   max_seurat_resolution = 5, ## change this to higher values
                                   FDRthresh = 0.05, # FDR threshold for statistical tests
                                   min_num_DE = 1,
                                   seurat_resolution = 0, # Starting resolution is this plus the jump value below.
                                   seurat_resolution_jump = 0.2,
                                   DE_bw_clust = TRUE)

### add the final clustering attributes to the meta data
sCVdata_list <- readRDS('Results/old_samples/sCVdata_list_old_data_final.rds')


########################


## creating the object
your_scRNAseq_data_object <- CreateSeuratObject(
  counts=assay_data,
  project = "MergedRats",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)

######### generating the embedding data frames
df_pca = Embeddings(merged_samples, 'pca')
rownames(df_pca) = corrected_UMIs

df_harmony = Embeddings(merged_samples, 'harmony')
rownames(df_harmony) = corrected_UMIs

df_umap = Embeddings(merged_samples, 'umap')
rownames(df_umap) = corrected_UMIs

df_tsne = Embeddings(merged_samples, 'tsne')
rownames(df_tsne) = corrected_UMIs


### adding dimension reduction embeddings
your_scRNAseq_data_object[["pca"]] <- CreateDimReducObject(embeddings = df_pca, 
                                                           key = "PCA_", 
                                                           assay = DefaultAssay(your_scRNAseq_data_object))

your_scRNAseq_data_object[["harmony"]] <- CreateDimReducObject(embeddings = df_harmony, 
                                                               key = "Harmony_", 
                                                               assay = DefaultAssay(your_scRNAseq_data_object))

your_scRNAseq_data_object[["umap"]] <- CreateDimReducObject(embeddings = df_umap, 
                                                            key = "UMAP_", 
                                                            assay = DefaultAssay(your_scRNAseq_data_object))

your_scRNAseq_data_object[["tsne"]] <- CreateDimReducObject(embeddings = df_tsne, 
                                                            key = "tSNE_", 
                                                            assay = DefaultAssay(your_scRNAseq_data_object))


your_scRNAseq_data_object$strain = merged_samples$strain
your_scRNAseq_data_object$orig.ident <- merged_samples$strain ## samples_name
your_scRNAseq_data_object$cluster = as.character(sCVdata_list$RNA_snn_res.1@Clusters)


#### new corrected set-2 objects
#saveRDS(your_scRNAseq_data_object, 'Objects/merged_samples_newSamples_MT-removed_labelCorrect.rds')
#saveRDS(your_scRNAseq_data_object, 'Results/new_samples/Immune_subclusters_labelCorrect.rds')
#saveRDS(your_scRNAseq_data_object, 'Results/new_samples/endothelial_subclusters_labelCorrect.rds')

## creating the meta.data dataframe
head(merged_samples@meta.data)
your_cluster_results <- data.frame(merged_samples@meta.data[,c(9:ncol(merged_samples@meta.data))]) # 7
rownames(your_cluster_results)  = rownames(merged_samples@meta.data)
head(your_cluster_results)

### calculating the differentially expressed marker genes
sCVdata_list <- CalcAllSCV(
  inD=your_scRNAseq_data_object,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.1, #gene filter - minimum detection rate
  testAll=F, #stop testing clusterings when no DE between clusters
  FDRthresh=0.05,
  calcSil=T, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T
)

'~/RatLiver/'
#### talk to Sonya about these scClusViz object ####
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included_labelCor.RData"
#_labelCor tag was not added although the labels have been refined:
new_data_scCLustViz_object_endothelial <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_EndothelialSub.RData"

new_data_scCLustViz_object <- "Results/new_samples/for_scClustViz.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_old_samples_highRes.RData"



for(i in 1:length(sCVdata_list)){
  names(sCVdata_list[[i]]@Clusters) = corrected_UMIs
}
save(your_scRNAseq_data_object,sCVdata_list,
      file=new_data_scCLustViz_object_endothelial) ## new data scClustViz object

your_scRNAseq_data_object$orig.ident <- merged_samples$orig.ident
save(your_scRNAseq_data_object,sCVdata_list,
      file=new_data_scCLustViz_object) ## old data scClustViz object


# This file can now be shared so anyone 
# can view your results with the Shiny app!
load(new_data_scCLustViz_object)
load(old_data_scClustViz_object)

load(new_data_scCLustViz_object_Immune)
load(new_data_scCLustViz_object_endothelial)


#### generating count table ####
df = data.frame(table(sCVdata_list$RNA_snn_res.1@Clusters))
df$Var1 = paste0('cluster_', as.character(df$Var1))
colnames(df) = c('set2-endothelial', '#cells')
gridExtra::grid.table(df)
####################


runShiny(
  ## write the path to the file bellow:
  filePath= new_data_scCLustViz_object,
  
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  
  annotationDB="org.Rn.eg.db",
  # This is an optional argument, but will add annotations.
  
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)


'~/RatLiver/Results/new_samples/rotated_loadings/rot_PC5_loadings.rnk'
'~/RatLiver/Results/new_samples/rotated_loadings_immune/rot_PC6_loadings.rnk'
merged_samples <- your_scRNAseq_data_object
#### check if removal of mt-genes would help with the annotation
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      clusters=paste0('', sCVdata_list$res.0.6@Clusters))

title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(values = colorPalatte)+
  theme(text = element_text(size=22),legend.title = element_blank())


df_tsne <- data.frame(tSNE_1=getEmb(merged_samples, 'tsne')[,1], 
                      tSNE_2=getEmb(merged_samples, 'tsne')[,2], 
                      clusters=paste0('cluster_', sCVdata_list$res.1@Clusters))

ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=clusters))+geom_point()+
  theme_classic()#+scale_color_manual(values = colorPalatte)

'~/RatLiver/Results/new_samples/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData'
'~/RatLiver/Results/new_samples/for_scClustViz_newSamples_MTremoved_EndothelialSub.RData'



