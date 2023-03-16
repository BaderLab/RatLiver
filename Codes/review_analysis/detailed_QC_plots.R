### loading the required libraries
source('Codes/Functions.R')
Initialize()


## Define cut-off values
MIT_CUT_OFF = 40
LIB_SIZE_CUT_OFF = 2000
NUM_GENES_DETECTED = 250
INPUT_NAME = 'rat_LEW_M09_WK_009_3pr_v3'# rat_DA_01_reseq 'rat_DA_M_10WK_003' # 'rat_DA_01_reseq' #rat_Lew_01 #'rat_Lew_02' 

##### TLH map
if(INPUT_NAME == 'rat_DA_01_reseq') {LIB_SIZE_CUT_OFF=1500; MIT_CUT_OFF=30}
if(INPUT_NAME == 'rat_DA_M_10WK_003') {LIB_SIZE_CUT_OFF=2000; MIT_CUT_OFF=20}

##### Immune enriched map
if(INPUT_NAME == 'rat_DA_M09_WK_008_3pr_v3') {LIB_SIZE_CUT_OFF=1000; MIT_CUT_OFF=50}
if(INPUT_NAME == 'rat_LEW_M09_WK_009_3pr_v3') {LIB_SIZE_CUT_OFF=1000; MIT_CUT_OFF=50}

# rat_LEW_M09_WK_009_3pr_v3: mito=40, library size= 1500


# rat_DA_01_reseq --> 50778
# rat_DA_M_10WK_003 --> 51232
# rat_Lew_01 --> 62362
# rat_Lew_02 --> 61898
# 50778 + 51232 + 62362 + 61898

#### new syncronized thresholds: 
MIT_CUT_OFF = 40
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250

TITLE = paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF)
#TITLE = paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', NUM_GENES_DETECTED)

OUTPUT_NAME = paste0('seur_QC_',INPUT_NAME,'_','mito_', MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.rds')
#OUTPUT_NAME = paste0('seur_QC_',INPUT_NAME,'_','mito_', MIT_CUT_OFF,'_numGene_',NUM_GENES_DETECTED,'.rds')


# ## import the data
input_from_10x <- paste0("Data/", INPUT_NAME,'/')
seur_raw <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                               min.cells=0,min.features=1, 
                               project = "snRNAseq")
dim(seur_raw)

seur_genes_df <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V2, col.name = 'symbol')
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V1, col.name = 'ensembl')

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, seur_raw[['RNA']]@meta.features$symbol )
seur_raw[['RNA']]@meta.features$symbol[mito_genes_index]
seur_raw[["mito_perc"]] <- PercentageFeatureSet(seur_raw, features = mito_genes_index)
summary(seur_raw[["mito_perc"]]$mito_perc )
libSize <- colSums(GetAssayData(seur_raw@assays$RNA))

##### calcularing the ribosamal fractionn as well

RIBO_PATTERN = '^Rp[sl]'
ribo_genes_index <- grep(pattern = RIBO_PATTERN, seur_raw[['RNA']]@meta.features$symbol )
seur_raw[['RNA']]@meta.features$symbol[ribo_genes_index]
seur_raw[["ribo_perc"]] <- PercentageFeatureSet(seur_raw, features = ribo_genes_index)

head(seur_raw@meta.data)



VlnPlot(seur_raw, features = c("nFeature_RNA", "nCount_RNA", "mito_perc", 'ribo_perc'), ncol = 4)
meta_qc_df = data.frame(seur_raw@meta.data, umi=INPUT_NAME)


to_drop_mito <- seur_raw$mito_perc > MIT_CUT_OFF
LIB_SIZE_CUT_OFF_MAX = 75000
to_drop_lib_size <- seur_raw$nCount_RNA < LIB_SIZE_CUT_OFF | seur_raw$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
to_drop_num_genes <- seur_raw$nFeature_RNA < NUM_GENES_DETECTED
df = data.frame(library_size= seur_raw$nCount_RNA, 
                mito_perc=seur_raw$mito_perc , 
                n_expressed=seur_raw$nFeature_RNA)

ggplot(data.frame(seur_raw$nCount_RNA), aes(seur_raw.nCount_RNA))+
  geom_histogram(bins = 60,color='black',fill='pink',alpha=0.5)+
  theme_bw()+ggtitle('library size for all cells (before filter)')+xlab('Library sizes')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(data.frame(seur_raw$nFeature_RNA), aes(seur_raw.nFeature_RNA))+
  geom_histogram(bins = 60,color='black',fill='blue',alpha=0.3)+
  theme_bw()+ggtitle('# expressed genes for all cells(before filtering)')+xlab('Number of expressed genes')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(data.frame(seur_raw$mito_perc), aes(seur_raw.mito_perc))+
  geom_histogram(bins = 60,color='black',fill='green',alpha=0.3)+
  theme_bw()+ggtitle('proportion of reads mapped to Mt genes(before filtering)')+xlab('Mitochondrial proportion (%)')+
  ylab('Number of cells')+labs(caption = INPUT_NAME)

ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  scale_color_viridis('#expressed\ngenes',direction = +1)+
  ggtitle(paste0('Mitochondrial transcript threshold: ', MIT_CUT_OFF,'\nLibrary size threshold: ', LIB_SIZE_CUT_OFF))


ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point(size=1.5,alpha=0.7)+scale_color_viridis()+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red", size=1)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=1)+
  #ggtitle('Library size and mitochondrial transcript cutoffs')+
  theme(axis.title = element_text(size = 17), axis.text = element_text(size = 16), plot.title=element_text(size=18))


ggplot(df, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('num detected genes')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  #ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', NUM_GENES_DETECTED, ' (before filter)'))+
  geom_vline(xintercept = NUM_GENES_DETECTED, linetype="dashed", color = "red3", size=0.5)
  


ggplot(df, aes(x=library_size, y=n_expressed, color=mito_perc))+geom_point()+labs(caption = INPUT_NAME)+
  theme_bw()+xlab('Library Size')+ylab('Number of expressed genes')+scale_color_viridis(option = 'magma')+ggtitle('before filter')

ggplot(df, aes(x=library_size, y=n_expressed))+geom_point(color='darkblue')+theme_bw()+xlab('Library Size')+
  ylab('Number of expressed genes')+geom_point(data=df[to_drop_mito,],pch=4,color="red")+labs(caption = INPUT_NAME)+
  scale_color_viridis(option='magma', direction = 1)+ggtitle('before filter, labels: high-mito cells')+labs(caption = INPUT_NAME)


seur <- seur_raw[,!to_drop_mito & !to_drop_lib_size]
show(seur)

df_filt = data.frame(library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA)
ggplot(df_filt, aes(x=library_size, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+labs(caption = INPUT_NAME)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))

ggplot(df_filt, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Number of detected genes')+ylab('Mitochondrial transcript percent')+labs(caption = INPUT_NAME)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', LIB_SIZE_CUT_OFF, ' (after filter)'))

