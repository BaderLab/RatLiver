source('Codes/Functions.R')
Initialize()
library(scDblFinder)
library(BiocParallel)
set.seed(123)

library('scDblFinder')

############
library(BiocParallel)
library(parallel)
#merged_data <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC.rds')
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
#merged_samples <- FindNeighbors(merged_samples,reduction="harmony",verbose=T)
Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)

sce = GetAssayData(merged_samples, 'counts')
sce = scDblFinder(sce, samples=merged_samples$sample_name, BPPARAM=MulticoreParam(detectCores()-2))


table(sce$scDblFinder.class)
doublet_df = data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1],
                        UMAP_2=getEmb(merged_samples, 'umap_h')[,2],
                        db=as.character(sce$scDblFinder.class),
                        umi=colnames(sce),
                        nuclear_fraction=merged_samples$nuclear_fraction,
                        cluster=merged_samples$SCT_snn_res.2.5,
                        Library_size=merged_samples$nCount_RNA,
                        mito_fraction=merged_samples$mito_perc)

ggplot2::ggplot(doublet_df, aes(UMAP_1, UMAP_2, color=db))+geom_point(alpha=0.5,size=1)+theme_classic()+
  scale_color_manual(values = c('red', 'grey86'))+
  theme(text = element_text(size=15),legend.title = element_blank())#

data.frame(table(doublet_df$db, as.character(merged_samples$SCT_snn_res.2.5)))

library(viridis)
ggplot(doublet_df, aes(UMAP_1, UMAP_2, color=nuclear_fraction))+geom_point(alpha=0.5,size=1.1)+theme_classic()+
  scale_color_viridis(option='magma',direction = -1)
ggplot(doublet_df, aes(UMAP_1, UMAP_2, color=mito_fraction))+geom_point(alpha=0.5,size=1.1)+theme_classic()+
  scale_color_viridis(option='viridis',direction = -1)
ggplot(doublet_df, aes(x=cluster, y=Library_size))+geom_boxplot()+theme_classic()




####################################################################################
##################### Running soupX using the raw merged and sub-file 
soupx_c = SoupChannel(tod=seur_merged, toc=seur_mergedSub, calcSoupProfile=T)


load("Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData")
seur = your_scRNAseq_data_object
immuneSub_df = data.frame(getEmb(seur, 'umap'))
immuneSub_df$strain = sapply(str_split(colnames(seur), '_'), '[[', 2)
### flipping the sample tag
immuneSub_df$strain = ifelse(immuneSub_df$strain == 'DA', 'LEW', 'DA')
umi_only = sapply(str_split(colnames(seur), '_'), '[[', 6)
umi_sample_info = ifelse(immuneSub_df$strain == "DA", 'rat_DA_M09_WK_008_', 'rat_LEW_M09_WK_009_')
#colnames(merged_samples) = paste0(umi_sample_info, umi_only)
immuneSub_df$umi = paste0(umi_sample_info, umi_only)


immuneSub_df = merge(immuneSub_df, doublet_df[,3:4], by.x='umi', by.y='umi', all.x=T)
head(immuneSub_df)
ggplot2::ggplot(immuneSub_df, aes(UMAP_1, UMAP_2, color=db))+geom_point(alpha=0.5,size=1)+theme_classic()+
  scale_color_manual(values = c('red', 'grey86'))+
  theme(text = element_text(size=15),legend.title = element_blank())#


####################################################################################
####################################################################################
####################################################################################
load("Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included_labelCor.RData")
seur = your_scRNAseq_data_object

immuneSub_df = data.frame(getEmb(seur, 'umap'))
immuneSub_df$cluster = as.character(sCVdata_list$res.1@Clusters)
immuneSub_df$umi = colnames(seur)

set1_info <- read.csv('figure_panel/set-2-immunesub-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info <- set1_info[1:14,]
set1_info$clusters = as.character(set1_info$clusters)
immuneSub_df = merge(immuneSub_df, set1_info[,1:3], by.x='cluster', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
immuneSub_df <- immuneSub_df[match(colnames(seur),immuneSub_df$umi),]
immuneSub_df$umi == colnames(seur)


library(UCell)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- length(names(table(seur$cluster))) #18 #14
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)


my.matrix <- GetAssayData(seur)
colnames(my.matrix) = immuneSub_df$umi
gene.sets <- list(Hep_signature = c('Gc','Serpina1','Ttr','Apoc3','Ifitm3','Apoa1','Fabp1','Ambp'))

scores <- ScoreSignatures_UCell(my.matrix, features=gene.sets)
immuneSub_df$Ucell_HepScore = scores[,1]
head(scores)
hist(scores[,1])


ggplot2::ggplot(immuneSub_df, aes(UMAP_1, UMAP_2, color=Ucell_HepScore))+
  geom_point(alpha=0.5,size=1)+theme_classic()+
  scale_color_viridis(direction = -1)+ggtitle('UCell Hep signature scores over\nset-2 immune-subset UMAP')
theme(text = element_text(size=15))#,legend.title = element_blank()

ggplot(immuneSub_df, aes(x=label, y=Ucell_HepScore, fill=label))+
  geom_boxplot()+theme_classic()+
  scale_fill_manual(values = mycolors)+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'))#,legend.title = element_blank()


HepMac_df = data.frame(umi=immuneSub_df$umi,
                       HepScore=immuneSub_df$Ucell_HepScore,
                       label=immuneSub_df$label)[as.character(immuneSub_df$label)=="Hep & non-Inf Mac (6)",]
head(HepMac_df)
ggplot(HepMac_df, aes(x=HepScore))+theme_classic()+
  geom_histogram(color="darkblue", fill="lightblue",aes(y = ..density..)) + 
  geom_density()+ggtitle('Ucell-Hep score distribution of the Mac&Hep(6) cluster\nthereshold selection')

thershold=0.75
HepMac_df$cells_2b_removed = HepMac_df$HepScore > thershold
UMIs_2bRemoved = HepMac_df$umi[HepMac_df$cells_2b_removed] #34 cells

seur$label = immuneSub_df$label
seur_nonHep = seur[,!colnames(seur)%in%UMIs_2bRemoved]
saveRDS(seur_nonHep, 'Results/new_samples/Immune_subclusters_labelCorrect_HepRemoved.rds')

