source('Codes/Functions.R')
Initialize()
library(remotes)
library(UCell)
library(ggplot2)
library(RColorBrewer)

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


PPARG_geneSets <- readRDS('~/PPARG_geneSets_converted2rat.rds')
PPARG_geneSets$`PPARG 20176806 ChIP-Seq MACROPHAGES Mouse`


my.matrix <- UCell::sample.matrix
my.matrix <- GetAssayData(your_scRNAseq_data_object)
gene.sets <- PPARG_geneSets
names(gene.sets) <- make.unique(sapply(strsplit(names(gene.sets), ' '), function(x) paste(x[1],x[4],sep = '-')))

### calculating enrichment score for each signature
scores <- data.frame(ScoreSignatures_UCell(my.matrix, features=gene.sets))

### percentage of genes included in the geneset which are present in the original dataset 
data.frame(length=sapply(gene.sets, function(a_geneset) length(a_geneset)),
           fraction=sapply(gene.sets, function(a_geneset) sum(a_geneset%in%rownames(my.matrix))/length(a_geneset)))


scores$strain <- sapply(strsplit(colnames(your_scRNAseq_data_object),'_'), function(x) x[2])
head(scores)

p1=ggplot(scores, aes(y=PPARG.3T3.L1_UCell, x=strain))+geom_violin(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)
p2=ggplot(scores, aes(y=PPARG.THIOMACROPHAGE_UCell, x=strain))+geom_violin(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)
p3=ggplot(scores, aes(y=PPARG.C3H10T1.2_UCell, x=strain))+geom_violin(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)
p4=ggplot(scores, aes(y=PPARG.3T3.L1.1_UCell, x=strain))+geom_violin(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)
p5=ggplot(scores, aes(y=PPARG.MACROPHAGES_UCell, x=strain))+geom_violin(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)
p6=ggplot(scores, aes(y=PPARG.3T3.L1.2_UCell, x=strain))+geom_violin(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)



######## selecting cells of interest
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rot_data <- readRDS('Results/new_samples/varimax_rotated_object_new.rds') ## MT-removed

rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
head(scores)
head(rotatedLoadings)

library(stats)
library(ggpubr)
kmeans_df = data.frame(scores$Varimax_1, scores$Varimax_15)
kmeans_res = kmeans(scale(kmeans_df),3)
kmeans_df$labels = as.character(kmeans_res$cluster)
ggplot(kmeans_df, aes(x=scores.Varimax_1,y=scores.Varimax_15,color=labels))+geom_point()+
  scale_color_brewer(palette = 'Set1')+ggtitle('Kmeans clustering based on VarPC-1 & 15')

kmeans_df <- readRDS('Results/kmeans_var15_1_oldsamples.rds')
high_score_cells = rownames(scores)[kmeans_df$labels %in% c('1','2')]



new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

### check the lables on the total umap
umap_df <- data.frame(Embeddings(your_scRNAseq_data_object, 'umap'),
                      Kmeans_label=as.character(kmeans_df$labels))

ggplot(umap_df, aes(UMAP_1, UMAP_2,color=Kmeans_label))+
  geom_point(alpha=0.6,size=2)+theme_classic()+
  scale_color_manual(values = c('orangered', 'deepskyblue2', 'forestgreen'))


my.matrix <- GetAssayData(your_scRNAseq_data_object)
my.matrix_sub = my.matrix[,colnames(my.matrix)%in%high_score_cells]
gmtFile = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'
geneSets <- getGmt(gmtFile)
names(geneSets)
query_genesets = names(geneSets)[grepl(pattern = 'lymphocyte migration', names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features=geneIds(geneSets[names(geneSets) %in% query_genesets])))
scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
colnames(scores_2)
p1=ggplot(scores_2, aes(y=POSITIVE.REGULATION.OF.LYMPHOCYTE.MIGRATION.GOBP.GO.2000403_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('POSITIVE.REGULATION.OF.LYMPHOCYTE.MIGRATION')+
  scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
p2=ggplot(scores_2, aes(y=LYMPHOCYTE.MIGRATION.GOBP.GO.0072676_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('LYMPHOCYTE.MIGRATION')+
  scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
p3=ggplot(scores_2, aes(y=REGULATION.OF.LYMPHOCYTE.MIGRATION.GOBP.GO.2000401_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('REGULATION.OF.LYMPHOCYTE.MIGRATION')+
  scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
gridExtra::grid.arrange(p1,p2,p3,nrow=1,ncol=3)



query_genesets = names(geneSets)[grepl(pattern = 'response to interferon-beta', names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features=geneIds(geneSets[names(geneSets) %in% query_genesets])))
colnames(scores_2) <- c('CELLULAR.RESPONSE.TO.INTERFERON.BETA', 'RESPONSE.TO.INTERFERON.BETA')
scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])

p1=ggplot(scores_2, aes(y=CELLULAR.RESPONSE.TO.INTERFERON.BETA, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('CELLULAR.RESPONSE.TO.INTERFERON.BETA')+
  scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
p2=ggplot(scores_2, aes(y=RESPONSE.TO.INTERFERON.BETA, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('RESPONSE.TO.INTERFERON.BETA')+
  scale_fill_brewer(palette = 'Dark2')+stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
gridExtra::grid.arrange(p1,p2,nrow=1,ncol=2)



geneSets <- readRDS('~/geneSets_converted2rat.rds') ### CHEA genesets converted to rat genes
a_gene_set_name <- 'SPI1' # 'HNF4A' # GATA1
query_genesets = names(geneSets)[grepl(pattern = a_gene_set_name, names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features= geneSets[names(geneSets) %in% query_genesets]))
colnames(scores_2) <- sapply(strsplit(colnames(scores_2),'\\.'),function(x) paste(x[1],x[5],x[6],sep = '_'))

scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
head(scores_2)

ggplot(scores_2, aes(y=SPI1_HL60_Human_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p1=ggplot(scores_2, aes(y=SPI1_K562_Human_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p2=ggplot(scores_2, aes(y=SPI1_THIOMACROPHAGE_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p3=ggplot(scores_2, aes(y=SPI1_GC_B, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
p4=ggplot(scores_2, aes(y=SPI1_ERYTHROLEUKEMIA_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
ggplot(scores_2, aes(y=SPI1_NB4_Human_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()
gridExtra::grid.arrange(p1,p2,p4,nrow=1,ncol=3)



a_gene_set_name <- 'HNF4A' # 'HNF4A' # GATA1
query_genesets = names(geneSets)[grepl(pattern = a_gene_set_name, names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, 
                                             features= geneSets[names(geneSets) %in% query_genesets]))
colnames(scores_2) <- sapply(strsplit(colnames(scores_2),'\\.'),function(x) paste(x[1],x[5],x[7],sep = '_'))

scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
head(scores_2)

ggplot(scores_2, aes(y=STAT4_TH1_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  scale_fill_brewer(palette = 'Dark2')+stat_compare_means()

