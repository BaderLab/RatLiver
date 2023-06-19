source('Codes/Functions.R')
Initialize()
library(remotes)
library(UCell)
library(ggplot2)
library(RColorBrewer)
library(stats)
library(ggpubr)



# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


var_num = 15
genes_num = 30

#### loading the data

## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
# merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)
merged_samples$clusters = merged_samples$SCT_snn_res.2.5




### loading the geneset
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rotatedLoadings <- rot_data$rotLoadings
head(rotatedLoadings)

### defining the signature based on varimax results
rotatedLoadings_ord = rotatedLoadings[order(rotatedLoadings[,var_num], decreasing = T),]
varimax_df = data.frame(genes=rownames(rotatedLoadings_ord) ,loading=rotatedLoadings_ord[,var_num])
gene_signature = varimax_df$genes[1:genes_num]

my.matrix <- GetAssayData(merged_samples)
gene.sets <- list(signature=gene_signature)

### calculating enrichment score for each signature
scores <- data.frame(ScoreSignatures_UCell(my.matrix, features=gene.sets))
head(scores)
sum(rownames(scores) != colnames(merged_samples))
scores.df = cbind(scores, merged_samples@meta.data, Embeddings(merged_samples, 'umap_h'))
head(scores.df)
ggplot(scores.df, aes(x=umap_h_1, y=umap_h_2, color=signature_UCell))+geom_point(alpha=0.4)+theme_classic()+
  scale_color_viridis('Strain Sig(PC15)\nenrichment score', option = 'inferno',direction = -1)+
  theme_classic()+
  xlab('UMAP_1')+ylab('UMAP_2')+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 13))+
  ggtitle('') 

table(scores.df$annot_TLH, scores.df$clusters)

ggplot(scores.df, aes(y=signature_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  scale_fill_brewer(palette = 'Dark2')+stat_compare_means()+theme_classic()+ggtitle('Complete maps')

ggplot(scores.df[scores.df$clusters==19,], aes(y=signature_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  scale_fill_brewer(palette = 'Dark2')+stat_compare_means()+theme_classic()+ggtitle('Macrophage clusters')+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 13))+ggtitle('')+ylab('Strain Sig (PC15) enrichment score')+xlab('')






