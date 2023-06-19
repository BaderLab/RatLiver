source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)
########################################################
############ Importing the nuc-seq data  #########

## importing the gene expression data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
# merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')

Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)
merged_samples$clusters = merged_samples$SCT_snn_res.2.5


cluster_num = '4'
gene_name = 'Cd5l'  #Sirpa
rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      clusters=merged_samples$clusters,
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      a_cluster=merged_samples$cluster==cluster_num,
                      umi=colnames(merged_samples))



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
  geom_point(alpha=0.6,size=1.2)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 

##### adding the final annotations to the dataframe
annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_May_31_2023.csv')
annot_info <- read.csv('figure_panel/Annotations_SingleNuc_Rat_June_8_2023.csv')
annot_info <- annot_info[1:34, 1:4]
colnames(annot_info)[1] = 'clusters'

################################################################
############### No need to run this for the May document
################################################################
#### annotation is not provided for the following clusters
clusters_not_included = names(table(merged_samples$SCT_snn_res.2.5))[!names(table(merged_samples$SCT_snn_res.2.5)) %in% 
                                               as.character(annot_info$clusters)]
  
clusters_not_included.df = data.frame(clusters=clusters_not_included,
                                      Annotation='Unknown',
                                      label='unknown',Markers=NA)
annot_info=rbind(annot_info, clusters_not_included.df)
#annot_info <- annot_info[1:18,] # 14
annot_info$clusters = as.character(annot_info$clusters)
################################################################


################################################################
########  adding color information to the data frame ##########
###############################################################

fc <- colorRampPalette(c("green", "darkgreen"))
annot_info.df = data.frame(table(annot_info$label))

fc <- colorRampPalette(c('pink1', 'palevioletred1'))
fc <- colorRampPalette(c('palevioletred1'))
cvHep_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Hep 1']) #cvHep

fc <- colorRampPalette(c('orchid1', 'orchid4'))
fc <- colorRampPalette(c('orchid3'))
Hep_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Hep 2'])

#fc <- colorRampPalette(c('palevioletred1', 'palevioletred1'))
fc <- colorRampPalette(c('magenta1', 'magenta3'))
fc <- colorRampPalette(c('maroon1'))
ppHep_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Hep 3']) #ppHep

fc <- colorRampPalette(c('grey70', 'grey20'))
fc <- colorRampPalette(c('mistyrose'))
unknown_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Unknown/High Mito'])

fc <- colorRampPalette(c('cyan')) #coral
chol_c = fc(annot_info.df$Freq[annot_info.df$Var1=='Cholangiocytes'])

infMac_c = '#117733'
nonInfMac_c = '#999933' #'#47A265'
LSEC_c = c('#FFE729','#FFAF13')
Stellate_c = '#6699CC'

color_df = data.frame(colors=c(cvHep_c, Hep_c, ppHep_c, chol_c, nonInfMac_c, infMac_c, Stellate_c, LSEC_c , unknown_c),
                      labels=c( rep('Hep 1',length(cvHep_c)), rep('Hep 2',length(Hep_c)), rep('Hep 3',length(ppHep_c)), 
                     rep('Cholangiocytes',length(chol_c)), 'Non-inf Macs', 'Inf Macs', 'Mesenchymal', 
                     rep('Endothelial',length(LSEC_c)), rep('Unknown/High Mito',length(unknown_c))))



color_df$labels == annot_info$label
show_col(color_df$colors)
annot_info2 = cbind(annot_info, color_df)

(0:33)[!(0:33) %in% as.numeric(annot_info$clusters)]
########## merging df_umap with annot_info data.frame
df_umap = merge(df_umap, annot_info2, by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
sum(df_umap$umi != colnames(merged_samples))
df_umap$label_clust = paste0(df_umap$label, ' (',df_umap$clusters, ')')


annot_info2$label_clust = paste0(annot_info2$label, ' (',annot_info2$clusters, ')')

Colors = annot_info2$colors
names(Colors) <- as.character(annot_info2$label)
names(Colors) <- as.character(annot_info2$label_clust)


#Colors = annot_info2$colors[match(levels(df_umap$clusters), annot_info2$clusters)]
#Colors = annot_info2$colors[match(levels(df_umap$clusters), annot_info2$clusters)]

library(scales)
show_col(Colors)

title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.8, size=1)+
  theme_classic()+scale_color_manual(values = Colors)+ #name='clusters',
  theme(text = element_text(size=15),legend.title = element_blank())#

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label_clust))+geom_point(alpha=0.8, size=1)+
  theme_classic()+scale_color_manual(name='label_clust',values = Colors)+
  theme(text = element_text(size=15),legend.title = element_blank())#

title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.7, size=2)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors)+
  theme(text = element_text(size=16),legend.title = element_blank())#

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.7, size=1.2)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors2)+
  theme(text = element_text(size=15),legend.title = element_blank())#

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=1, size=3)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors2)+
  theme(text = element_text(size=15),legend.title = element_blank())#




library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
################################################################
################ Generating a contribution plot ################
################################################################

df <- data.frame(sample_type = merged_samples$sample_name, 
                 cluster = as.character(df_umap$label))

#### based on sample-type ##### 
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")
counts$cluster = gsub(x =counts$cluster , 'Inflammatory', 'Inf')
counts$cluster = gsub(x =counts$cluster , 'inflammatory', 'inf')


counts$cluster= factor(counts$cluster, levels = as.character(Colors) ) 

counts$strain=sapply(strsplit(counts$sample_type,'-') ,'[[', 1)
#write.csv(counts, 'figure_panel/set1-counts.csv')


ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  #scale_fill_manual(values=Colors)+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  #scale_x_discrete(labels = xlabel2)+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )

ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Fraction of sample per cell type (%)')+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  #scale_fill_brewer(palette = "Set3")+
  #scale_fill_manual(values=mycolors)+
  #scale_x_discrete(labels = immune_sub_labels) #scales::wrap_format(9)
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12.5,angle=90,color='black'),
        legend.title = element_blank()) +  
  xlab('')





