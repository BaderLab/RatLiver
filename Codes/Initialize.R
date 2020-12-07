## Run this script as: 
# Rscript Codes/Init.R 'rat_Rnor' 50 1500

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

INPUT_NAME = args[1] 
MIT_CUT_OFF = as.numeric(args[2])
LIB_SIZE_CUT_OFF = as.numeric(args[3])
PC_NUMBER = 18
MADS_CUT_OFF = 12


source('Codes/Functions.R')
Initialize()

## Define cut-ff values
MIT_CUT_OFF = 50
LIB_SIZE_CUT_OFF = 1000
INPUT_NAME = 'rat_DA_M09_WK_008_3pr_v3' #'rat_LEW_M09_WK_009_3pr_v3'
#'rat_DA_M09_WK_008_3pr_v3'# 'rat_DA_01_reseq' 'rat_DA_M_10WK_003' 'rat_LEW_M09_WK_009_3pr_v3'

# MIT_CUT_OFF = median(seur$mito_perc) + mad(seur$mito_perc) * MADS_CUT_OFF
TITLE = paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF)
OBJ_NAME_NORM = paste0('1.seur_normed_',INPUT_NAME,'_','mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.rds')
OBJ_NAME_DIM_RED = paste0('2.seur_dimRed_',INPUT_NAME,'_','mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.rds')

OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))

# ## import the data
input_from_10x <- paste0("Data/", INPUT_NAME,'/')
# 
# ##### alternative input files
# emptyDrops_output <- readRDS('SoupX/rat_DA_01_emptyDrops_table.rds')
# seur <- CreateSeuratObject(emptyDrops_output)
# 
# SoupX_output <- readRDS('output_soupX_estimatedCont_10Xfilter.rds')
# seur <- CreateSeuratObject(SoupX_output)

########################
seur <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 1),
                           min.cells=0,min.features=1, 
                           project = "snRNAseq")

dim(seur)
seur_genes_df <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)
seur[['RNA']] <- AddMetaData(seur[['RNA']], seur_genes_df$V2, col.name = 'symbol')

MIT_PATTERN = '^Mt-'
if(INPUT_NAME == 'mouse') {MIT_PATTERN = '^mt-'} 
mito_genes_index <- grep(pattern = MIT_PATTERN, seur[['RNA']]@meta.features$symbol )
# mito_genes_index <- grep(pattern = MIT_PATTERN, rownames(seur) )

# pig_mit_genes_ensembl <- readRDS('Data/pig_mit_genes_ensembl.rds')
# mito_genes_index <- which(rownames(seur) %in% pig_mit_genes_ensembl)
seur[['RNA']]@meta.features$symbol[mito_genes_index]
seur[["mito_perc"]] <- PercentageFeatureSet(seur, features = mito_genes_index)
summary(seur[["mito_perc"]]$mito_perc )

## adding ensemble id as a meta data to the object

getHead(GetAssayData(seur@assays$RNA))
libSize <- colSums(GetAssayData(seur@assays$RNA))
summary(libSize)


print(paste0('Total number of cells: ', ncol(seur)))

to_drop_mito <- seur$mito_perc > MIT_CUT_OFF
print(paste0('to_drop_mito: ',sum(to_drop_mito)))
print(paste0('to_drop_mito percentage: ', round(sum(to_drop_mito)*100/ncol(seur),2) ))

LIB_SIZE_CUT_OFF_MAX = 75000
to_drop_lib_size <- seur$nCount_RNA < LIB_SIZE_CUT_OFF | seur$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
print(paste0('to_drop_lib_size: ', sum(to_drop_lib_size)))
print(paste0('to_drop_lib_size percentage: ', round( sum(to_drop_lib_size)*100/ncol(seur),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_lib_size & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_lib_size & !to_drop_mito)*100/ncol(seur),2) ))

df = data.frame(library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA)



## Visualization of QC metrics
pdf(paste0('Results/',INPUT_NAME,'/QC/QC_',INPUT_NAME,'_',
           'mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.pdf'))

ggplot(data.frame(seur$nCount_RNA), aes(seur.nCount_RNA))+
  geom_histogram(bins = 60,color='black',fill='pink',alpha=0.5)+
  theme_bw()+ggtitle('library size for all cells (before filter)')+xlab('Library sizes')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(data.frame(seur$nFeature_RNA), aes(seur.nFeature_RNA))+
  geom_histogram(bins = 60,color='black',fill='blue',alpha=0.3)+
  theme_bw()+ggtitle('# expressed genes for all cells(before filtering)')+xlab('Number of expressed genes')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(data.frame(seur$mito_perc), aes(seur.mito_perc))+
  geom_histogram(bins = 60,color='black',fill='green',alpha=0.3)+
  theme_bw()+ggtitle('proportion of reads mapped to Mt genes(before filtering)')+xlab('Mitochondrial proportion (%)')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(df, aes(x=library_size, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (before filter)'))

ggplot(df, aes(x=library_size, y=n_expressed, color=mito_perc))+geom_point()+labs(caption = INPUT_NAME)+
  theme_bw()+xlab('Library Size')+ylab('Number of expressed genes')+scale_color_viridis(option = 'magma')+ggtitle('before filter')

ggplot(df, aes(x=library_size, y=n_expressed))+geom_point(color='darkblue')+theme_bw()+xlab('Library Size')+
  ylab('Number of expressed genes')+geom_point(data=df[to_drop_mito,],pch=4,color="red")+labs(caption = INPUT_NAME)+
  scale_color_viridis(option='magma', direction = 1)+ggtitle('before filter, labels: high-mito cells')+labs(caption = INPUT_NAME)


seur <- seur[,!to_drop_mito & !to_drop_lib_size]
show(seur)
dim(seur)

df_filt = data.frame(library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA)
ggplot(df_filt, aes(x=library_size, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+labs(caption = INPUT_NAME)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))



## Normalization
### SCTransform

# seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
# dim(seur[['RNA']]@data)
# ## Finding variable genes
# seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
# head(seur[['RNA']]@var.features)
# ## scaling data
# seur <- ScaleData(seur, features = rownames(seur))

## alternative: 
seur <- SCTransform(seur,conserve.memory=F,verbose=T,return.only.var.genes=F,
                    variable.features.n = nrow(seur[['RNA']]))
dim(seur)
# saveRDS(seur, paste0('objects/',INPUT_NAME,'/',OBJ_NAME_NORM))

##  PCA
seur <- RunPCA(seur,verbose=T)
plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

dev.off()

PC_NUMBER = 22

pdf(paste0('Results/',INPUT_NAME,'/QC/plots_',INPUT_NAME,'_','mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.pdf'))

## tSNE
seur <- RunTSNE(seur,dims=1:PC_NUMBER,reduction="pca",perplexity=30)

TITLE_tsne = paste0('tSNE (',TITLE,')')
df_tsne <- data.frame(tSNE_1=getEmb(seur, 'tsne')[,1], tSNE_2=getEmb(seur, 'tsne')[,2], 
                      library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA) 

ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=library_size))+geom_point()+theme_bw()+scale_color_viridis(direction= -1)+ggtitle(TITLE_tsne)+labs(caption = INPUT_NAME)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=mito_perc))+geom_point()+theme_bw()+scale_color_viridis(direction= -1)+ggtitle(TITLE_tsne)+labs(caption = INPUT_NAME)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=n_expressed))+geom_point()+theme_bw()+scale_color_viridis(direction= -1)+ggtitle(TITLE_tsne)+labs(caption = INPUT_NAME)


##UMAP
seur <- RunUMAP(seur,dims=1:PC_NUMBER, reduction="pca")

saveRDS(seur, paste0('objects/',INPUT_NAME,'/',OBJ_NAME_DIM_RED))

TITLE_umap = paste0('UMAP (',TITLE,')')
df_umap <- data.frame(UMAP_1=getEmb(seur, 'umap')[,1], UMAP_2=getEmb(seur, 'umap')[,2], 
                      library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA )

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)


dev.off()

# objects/rat_LEW_M09_WK_009_3pr_v3/2.seur_dimRed_rat_LEW_M09_WK_009_3pr_v3_mito_50_lib_1000.rds
# objects/rat_DA_M09_WK_008_3pr_v3/2.seur_dimRed_rat_DA_M09_WK_008_3pr_v3_mito_50_lib_1000.rds

