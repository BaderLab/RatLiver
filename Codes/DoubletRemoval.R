source('Codes/Functions.R')
Initialize()
library(scDblFinder)
library(BiocParallel)
set.seed(123)

BiocManager::install('scDblFinder')

############
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

load(new_data_scCLustViz_object)
load(old_data_scClustViz_object)

merged_samples = your_scRNAseq_data_object

merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))

selected_UMIs = colnames(merged_samples)

######### importing the new samples raw data
seur_DA <- CreateSeuratObject(counts=Read10X('Data/rat_DA_M09_WK_008_3pr_v3/', gene.column = 2),min.cells=0,min.features=1, project = "snRNAseq")
seur_LEW <- CreateSeuratObject(counts=Read10X('Data/rat_LEW_M09_WK_009_3pr_v3/', gene.column = 2),min.cells=0,min.features=1, project = "snRNAseq")

seur_merged = merge(seur_DA, seur_LEW, add.cell.ids = c('rat_DA_M09_WK_008', 'rat_LEW_M09_WK_009'), 
      project = "rat_data", merge.data = TRUE)

###########################
######### importing the old samples raw data

seur_DA_1 <- CreateSeuratObject(counts=Read10X('Data/rat_DA_01_reseq/', gene.column = 2), min.cells=0,min.features=1, project = "snRNAseq")
seur_DA_2 = CreateSeuratObject(counts=Read10X('Data/rat_DA_M_10WK_003/', gene.column = 2), min.cells=0,min.features=1, project = "snRNAseq")
seur_LEW_1 <- CreateSeuratObject(counts=Read10X('Data/rat_Lew_01/', gene.column = 2), min.cells=0,min.features=1, project = "snRNAseq")
seur_LEW_2 <- CreateSeuratObject(counts=Read10X('Data/rat_Lew_02/', gene.column = 2), min.cells=0,min.features=1, project = "snRNAseq")

seur_merged = merge(seur_DA_1, c(seur_DA_2, seur_LEW_1, seur_LEW_2), 
                    add.cell.ids = c('rat_DA_01_reseq', 'rat_DA_M_10WK_003', 'rat_Lew_01', 'rat_Lew_02'), 
                    project = "rat_data", merge.data = TRUE)
###########################

#### all the cell names in the merged_samples (after filter) need to be included in the seur_merged object
sum(!colnames(merged_samples) %in% colnames(seur_merged))
seur_mergedSub = seur_merged[,colnames(seur_merged) %in% selected_UMIs]
dim(seur_mergedSub)

sce = GetAssayData(seur_mergedSub)
sce <- scDblFinder(sce, 
                   samples=sapply(str_split(colnames(seur_mergedSub), '_'), function(x) paste(x[-length(x)],collapse = '_')), 
                   BPPARAM=MulticoreParam(detectCores()-2))

table(sce$scDblFinder.class)
doublet_df = data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1],
                        UMAP_2=getEmb(merged_samples, 'umap')[,2],
                        db=as.character(sce$scDblFinder.class),
                        umi=colnames(sce))
ggplot2::ggplot(doublet_df, aes(UMAP_1, UMAP_2, color=db))+geom_point(alpha=0.5,size=1)+theme_classic()+
  scale_color_manual(values = c('red', 'grey86'))+
  theme(text = element_text(size=15),legend.title = element_blank())#

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

