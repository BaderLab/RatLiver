# In this script, we will read the 10X data files into session and
# we will filter the data based on the fraction of mitochondrial transcripts and 
# also library size, then we will scale and normalize the data and
# perform dimension reduction. Finally, we'll visualize the data using tsne and umap

## there are 6 datasets included in this analysis:
# rat_LEW_M09_WK_009_3pr_v3, rat_DA_M09_WK_008_3pr_v3, 
# rat_DA_M_10WK_003, rat_DA_01_reseq , rat_Lew_01 , rat_Lew_02

# rat_DA_01_reseq: mito=30, library size= 1500
# rat_DA_M_10WK_003: mito=20, library size=2000
# rat_DA_M09_WK_008_3pr_v3: mito=50,  library size=1000

# rat_Lew_01: mito=40, library size= 2000
# rat_Lew_02: mito=40, library size= 2000
# rat_LEW_M09_WK_009_3pr_v3: mito=40, library size= 1500


# Tallulah: 50% mt threshold across all my samples plus >=250 detected genes

### loading the required libraries
source('Codes/Functions.R')
Initialize()


## Define cut-off values
MIT_CUT_OFF = 40
LIB_SIZE_CUT_OFF = 2000
NUM_GENES_DETECTED = 250
INPUT_NAME = 'rat_DA_01_reseq'# 'rat_DA_M_10WK_003' # 'rat_DA_01_reseq' #'rat_Lew_02' 

if(INPUT_NAME == 'rat_DA_01_reseq') {LIB_SIZE_CUT_OFF=1500; MIT_CUT_OFF=30}
if(INPUT_NAME == 'rat_DA_M_10WK_003') {LIB_SIZE_CUT_OFF=2000; MIT_CUT_OFF=20}
if(INPUT_NAME == 'rat_DA_M09_WK_008_3pr_v3') {LIB_SIZE_CUT_OFF=1000; MIT_CUT_OFF=50}
if(INPUT_NAME == 'rat_LEW_M09_WK_009_3pr_v3') {LIB_SIZE_CUT_OFF=1500; MIT_CUT_OFF=40}


#### new syncronized thresholds: 
MIT_CUT_OFF = 40
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250

TITLE = paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF)
#TITLE = paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', NUM_GENES_DETECTED)

OUTPUT_NAME = paste0('seur_QC_',INPUT_NAME,'_','mito_', MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.rds')
#OUTPUT_NAME = paste0('seur_QC_',INPUT_NAME,'_','mito_', MIT_CUT_OFF,'_numGene_',NUM_GENES_DETECTED,'.rds')

OUTPUT_NAME

# ## import the data
input_from_10x <- paste0("Data/", INPUT_NAME,'/')
seur_raw <- CreateSeuratObject(counts=Read10X(input_from_10x, gene.column = 2),
                           min.cells=0,min.features=1, 
                           project = "snRNAseq")

seur_genes_df <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V2, col.name = 'symbol')
seur_raw[['RNA']] <- AddMetaData(seur_raw[['RNA']], seur_genes_df$V1, col.name = 'ensembl')

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, seur_raw[['RNA']]@meta.features$symbol )
seur_raw[['RNA']]@meta.features$symbol[mito_genes_index]
seur_raw[["mito_perc"]] <- PercentageFeatureSet(seur_raw, features = mito_genes_index)
summary(seur_raw[["mito_perc"]]$mito_perc )
libSize <- colSums(GetAssayData(seur_raw@assays$RNA))


print(paste0('Total number of cells: ', ncol(seur_raw)))

to_drop_mito <- seur_raw$mito_perc > MIT_CUT_OFF

print(paste0('to_drop_mito: ',sum(to_drop_mito)))
print(paste0('to_drop_mito percentage: ', round(sum(to_drop_mito)*100/ncol(seur_raw),2) ))

LIB_SIZE_CUT_OFF_MAX = 75000
to_drop_lib_size <- seur_raw$nCount_RNA < LIB_SIZE_CUT_OFF | seur_raw$nCount_RNA > LIB_SIZE_CUT_OFF_MAX
print(paste0('to_drop_lib_size: ', sum(to_drop_lib_size)))
print(paste0('to_drop_lib_size percentage: ', round( sum(to_drop_lib_size)*100/ncol(seur_raw),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_lib_size & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_lib_size & !to_drop_mito)*100/ncol(seur_raw),2) ))



to_drop_num_genes <- seur_raw$nFeature_RNA < NUM_GENES_DETECTED
print(paste0('to_drop_num_genes: ', sum(to_drop_num_genes)))
print(paste0('to_drop_num_genes percentage: ', round( sum(to_drop_num_genes)*100/ncol(seur_raw),2) ))

print(paste0('remained after both filters: ', sum(!to_drop_num_genes & !to_drop_mito)))
print(paste0('Percentage of remained after both filters: ', 
             round(  sum(!to_drop_num_genes & !to_drop_mito)*100/ncol(seur_raw),2) ))


df = data.frame(library_size= seur_raw$nCount_RNA, 
                mito_perc=seur_raw$mito_perc , 
                n_expressed=seur_raw$nFeature_RNA)


dir.create(paste0('Plots/QC/',INPUT_NAME))
## Visualization of QC metrics
pdf(paste0('Plots/QC/',INPUT_NAME,'/QC_',INPUT_NAME,'_',
           'mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.pdf'))

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

ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF, ' (before filter)'))


ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point(size=,alpha=0.7)+scale_color_viridis()+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  ggtitle('Library size and mitochondrial transcript cutoffs')



ggplot(df, aes(x=n_expressed, y=mito_perc, color=library_size))+geom_point()+scale_color_viridis()+
  theme_bw()+xlab('num detected genes')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+labs(caption = INPUT_NAME)+
  geom_vline(xintercept = NUM_GENES_DETECTED, linetype="dashed", color = "red3", size=0.5)+
  ggtitle(paste0('mito threshold: ', MIT_CUT_OFF,' , num detected genes threshold: ', NUM_GENES_DETECTED, ' (before filter)'))


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



## Normalization
### SCTransform
seur <- SCTransform(seur,conserve.memory=F,verbose=T,return.only.var.genes=F,
                    variable.features.n = nrow(seur[['RNA']]), 
                    ## according to the paper scaling is not recommended prior to PCA:
                    ## https://www.biorxiv.org/content/10.1101/576827v2
                    do.scale = FALSE, ### default value 
                    do.center = TRUE) ### default value 

dim(seur)

saveRDS(seur, paste0('Objects/', INPUT_NAME, '/',OUTPUT_NAME))





##  PCA
seur <- RunPCA(seur,verbose=T)
plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

dev.off()

PC_NUMBER = 18


pdf(paste0('Plots/QC/',INPUT_NAME,'/',INPUT_NAME,'_','mito_',MIT_CUT_OFF,'_lib_',LIB_SIZE_CUT_OFF,'.pdf'))
#pdf(paste0('Plots/QC/',INPUT_NAME,'/',INPUT_NAME,'_','mito_',MIT_CUT_OFF,'_numGenes_',NUM_GENES_DETECTED,'.pdf'))

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

dir.create(paste0('Objects/',INPUT_NAME))
saveRDS(seur, paste0('Objects/',INPUT_NAME,'/',OUTPUT_NAME))
seur <- readRDS(paste0('Objects/',INPUT_NAME,'/',OUTPUT_NAME))


TITLE_umap = paste0('UMAP (',TITLE,')')
df_umap <- data.frame(UMAP_1=getEmb(seur, 'umap')[,1], 
                      UMAP_2=getEmb(seur, 'umap')[,2], 
                      library_size= seur$nCount_RNA, 
                      mito_perc=seur$mito_perc , 
                      n_expressed=seur$nFeature_RNA )

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)

dev.off()


data <- readRDS('Objects/rat_DA_M09_WK_008_3pr_v3/seur_QC_rat_DA_M09_WK_008_3pr_v3_mito_50_lib_1000.rds')
dim(data)
summary(data$mito_perc)
summary(data$nCount_RNA)
summary(data$nFeature_RNA)


