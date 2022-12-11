source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()
library(plyr)
library(stats)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)


### strain gene figures are in the consistent_strain_gene.R
#### varimax-cluster correlations in the varimax_cluster_cor.R script

############
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

#new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
new_data_scCLustViz_object_Immune <- "~/RatLiver/Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included_labelCor.RData"
new_data_scCLustViz_object_endothelial <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_EndothelialSub.RData"

#########
load(new_data_scCLustViz_object)
load(old_data_scClustViz_object)

load(new_data_scCLustViz_object_Immune)
load(new_data_scCLustViz_object_endothelial)

#### the new immune subcluster data ####
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)
merged_samples$cluster = as.character(sCVdata_list$res.1@Clusters)
merged_samples$cluster = as.character(sCVdata_list$res.1@Clusters)
merged_samples$cluster = as.character(sCVdata_list$RNA_snn_res.1@Clusters)

table(merged_samples$orig.ident)
merged_samples$sample_name = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'LEW-1', 'LEW-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)

merged_samples$EndoSub = ifelse(colnames(merged_samples) %in% 
                                    colnames(your_scRNAseq_data_object), 'Endothelial\ncells', 'other')


cell_info_df = data.frame(merged_samples@meta.data, cell_id=colnames(merged_samples)) 
saveRDS(cell_info_df, 'Results/old_samples/MT_info/cell_info_df.rds')

######################################
############## Checking QC variables 
######################################
table(merged_samples$sample_name)
sample_name = 'LEW-1'
a_sample = merged_samples[,merged_samples$sample_name==sample_name]

ncol(a_sample)
### N genes per sample
sum(rowSums(a_sample) > 0 )
#### mmedian number of expressed genes per cell
median(colSums(GetAssayData(a_sample)>0))

#############################################################################
####### Evaluating the number of Ptprc+ cells in set-1 map vs set2 ##########
############################################################ 
##########  new immune enriched map : num 914, percentage:  0.24 of cells were Ptprc+
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0000  0.0000  0.3258  0.0000  2.7726

##########  old map : num 969, percentage:  0.04 of cells were Ptprc+
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00000 0.00000 0.00000 0.03149 0.00000 1.60944 
############################################################
Ptprc_exp = GetAssayData(merged_samples)['Ptprc',]
round(sum(Ptprc_exp>0)/length(Ptprc_exp),2) 
summary(Ptprc_exp)
summary(Ptprc_exp[Ptprc_exp>0])

##################################################################

RT_genes = rownames(merged_samples)[grep('Ptgdr', rownames(merged_samples))]


cluster_num = '4'
gene_name = 'Sirpa'  #Sirpa
rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      a_cluster=merged_samples$cluster==cluster_num)


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
  geom_point(alpha=0.6,size=1.2)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_cluster))+
  geom_point(alpha=0.6,size=1.2)
#########################################################

                      
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      label=merged_samples$EndoSub,
                      Ptprc=GetAssayData(merged_samples)['Ptprc',],
                      Lyve1=GetAssayData(merged_samples)['Lyve1',],
                      Eng=GetAssayData(merged_samples)['Eng',])
                      
#### check if removal of mt-genes would help with the annotation
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      #clusters=paste0('', sCVdata_list$res.0.6@Clusters),
                      clusters=merged_samples$cluster,
                      sample=merged_samples$sample_name,
                      strain=merged_samples$strain,
                      Ptprc=GetAssayData(merged_samples)['Ptprc',], 
                      umi=colnames(merged_samples), 
                      libSize=merged_samples$nFeature_RNA)

df_umap$Cd14 <-  GetAssayData(merged_samples)['Cd14',]
##### adding the final annotations to the dataframe
set1_info <- read.csv('figure_panel/set-1-final-info.csv')
set1_info <- read.csv('figure_panel/set-2-final-info-updated.csv')
set1_info <- read.csv('figure_panel/set-2-Endosub-final-info.csv')
set1_info <- read.csv('figure_panel/set-2-immunesub-final-info-updated.csv')
colnames(set1_info)[1] = 'clusters'
set1_info$Annotation = gsub(pattern = 'Naive', replacement ='ab' ,set1_info$Annotation)
#set1_info <- set1_info[1:18,] # 14
set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T, order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
sum(df_umap$umi != colnames(merged_samples))

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

# Define the number of colors you want
nb.cols <- length(names(table(merged_samples$cluster))) #18 #14
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
show_col(c(mycolors,colorBlindGrey8,safe_colorblind_palette))
matrix(c(mycolors,colorBlindGrey8,safe_colorblind_palette),nrow=6,ncol=7)

### Pink, brown
show_col(c('#F781BF','#D36E7C' , '#9E5198','#AA4499',   '#984560','#5E3C99','pink2','#882255','#661100'))
### yellow
show_col(c('#FFE729', '#FFAF13',  '#E66101', '#F87B0A'))
### green
show_col(c('#117733', '#44AA99', '#999933', '#47A265', '#3D8D95'))
### blue
show_col(c('#6699CC', '#5684E9', '#0072B2', '#332288'))
## red
c('#E41A1C')

show_col(c('#999999', 'black'))
#### Set1 map
mycolors2 = c(c('#F781BF','#D36E7C' , '#9E5198','#882255','#661100','#AA4499',   '#984560','pink2','#5E3C99'),
              '#117733' ,
              '#FFE729','#FFAF13', 
              '#E41A1C', 
              '#999933', '#47A265',
              '#6699CC', '#0072B2')

mycolors2 = c('#999999',
              '#E41A1C',
              '#F781BF', '#D36E7C', '#9F5198','#661100', '#984560','pink2','#5E3C99',
              '#117733',
              '#FFE729','#FFAF13',
              "#F87B0A",#mature B cell,
              "#D55E00", #'#E66101',#naive T cell,
              '#999933', '#47A265',
              '#3D8D95', # pDC 
              "#5684E9")


title = '' # 'Mt-genes removed'
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point(alpha=0.8, size=2)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors2)+
  theme(text = element_text(size=15),legend.title = element_blank())#

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point(alpha=0.8, size=1)+
  theme_classic()+scale_color_manual(name='clusters',values = mycolors)+
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

df_umap$label_hep <- ifelse(substr(df_umap$label, 1, 3)=='Hep', df_umap$label, 'Other')
df_umap$label_lsec <- ifelse(substr(df_umap$label, 1, 4)=='LSEC', df_umap$label, 'Other')
df_umap$label_stellate <- ifelse(substr(df_umap$label, 1, 8)=='Stellate', df_umap$label, 'Other')
df_umap$label_nonInfMac <- ifelse(substr(df_umap$label, 1, 7)=='Non-Inf', df_umap$label, 'Other')


#### calculating the percentage of Heps 
num_hep_clust = length(table(df_umap$label_hep))-1
round(sum(table(df_umap$label_hep)[1:num_hep_clust])/sum(table(df_umap$label_hep)),4)* 100


names(table(df_umap$label))
head(df_umap)
df_umap$num_genes = merged_samples$nFeature_RNA
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=num_genes))+geom_point(alpha=0.7, size=1)+
  theme_classic()+scale_color_viridis('#expressed\ngenes',direction = +1)+
  theme(text = element_text(size=17))#



nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(8, "Pastel1"))(nb.cols) #Pastel1

mycolors = c( "#FF9999", "#CC79A7","#56B4E9","#0072B2")


show_col(c(mycolors, "#D55E00"))
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample))+geom_point(alpha=1, size=3)+
  theme_classic()+
  #scale_color_manual(values = c("#E69F00", "#009E73"))+
  #scale_color_manual(values = c("pink1", "skyblue1"))+
  #scale_color_brewer(palette = "Set3")+
  scale_color_manual(values=mycolors)+
  theme(text = element_text(size=19),legend.title = element_blank())

##### Marker expression over UMAP

Stellate_markers = c('Calcrl', 'Ecm1', 'Igfbp7', 'Rbp1', 'Co3a1', 'Bgn', 
                     'Igfbp3',"Colec11", "Colec10", "Colec12")

Hep_markers = c('Cyp7a1', 'Oat', 'Scd', 'Axin2', 'Glul', 'Hal', 'Orm1')

pdf('Plots/Set1-stellateMarkers.pdf', height = 7, width = 8)
pdf('Plots/Set1-HepMarkers.pdf', height = 7, width = 8)

for(i in 1:length(Hep_markers)){
  marker_name = Hep_markers[i] # Cd11b/Itgam not included
  if(!marker_name %in% rownames(merged_samples)) next
  df_umap$marker=GetAssayData(merged_samples)[marker_name,]
  p=ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=marker))+geom_point(alpha=0.4,size=0.6)+
    scale_color_viridis(direction = -1)+ggtitle(marker_name)+theme_classic()+
    theme(text = element_text(size=17.5),plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
  print(p)
}
dev.off()


##### checking one specific cluster
a_cluster = as.character(c(7, 8, 10, 11, 12, 17)) 
a_cluster = as.character(c(1, 2, 3, 9, 14, 15)) 

df_umap$a_cluster = ifelse(df_umap$clusters %in% a_cluster,  df_umap$clusters, 'other')

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_cluster))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(name='Cluster',values = c('black', 'grey'))+
  theme(text = element_text(size=22))+#,legend.title = element_blank()
  ggtitle(paste0('Cluster-', a_cluster, ' cells over Set-2 UMAP'))


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_cluster))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(name='Clusters',values = c(rep('black', 6), 'grey'))+
  theme(text = element_text(size=22))+#,legend.title = element_blank()
  ggtitle(paste0('Cluster 1, 2, 3, 9, 14, 15 cells\nover Set-2 UMAP'))

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.6, size=4.5)+
  theme_classic()+scale_color_manual(values = c('black', 'grey'))+
  theme(text = element_text(size=22),legend.title = element_blank())+
  ggtitle(paste0('Endothelial subclusters\nSet-2 UMAP'))


################
df <- data.frame(sample_type = merged_samples$sample_name, 
                 cluster = as.character(df_umap$label))
rownames(df) = NULL
head(df)

nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(8, "Pastel1"))(nb.cols)

####### set1 labels
xlabel = c("Hep (0)","Hep (1)","Hep (2)","LSEC (3)","Hep (4)","Non-Inf Mac (5)",
           "Hep (6)","Stellate (7)", "Hep (8)","Inf Mac (9)", "Non-Inf Mac (10)", 
           "LSEC (11)", "Hep (12)", "NK-like & T cells (13)", "Stellate (14)", "Hep (15)", "Hep (16)" )  

set_1_cl_ord = c('Hep (2)'	, 'Hep (4)',  'Hep (8)',  'Hep (16)',  'Hep (1)',  
                 'Hep (6)',  'Hep (0)',  'Hep (12)','Hep (15)',	'LSEC (3)', 
                 'LSEC (11)',	'Stellate (7)','Stellate (14)',	'Non-Inflammatory Mac (5)',		
                 'Non-Inflammatory Mac (10)', 'Inflammatory Mac (9)',	'NK-like and T cells (13)')

xlabel2 = c('Hep (2)'	, 'Hep (4)',  'Hep (8)',  'Hep (16)',  'Hep (1)',  
                 'Hep (6)',  'Hep (0)',  'Hep (12)','Hep (15)',	'LSEC (3)', 
                 'LSEC (11)',	'Stellate (7)','Stellate (14)',	'NonInf Mac(5)',		
                 'NonInf Mac(10)', 'Inf Mac (9)',	'NK-like &\nT cells (13)')

##### set-2 immune-sub labels

immune_sub_labels = levels(counts_norm$cluster)
immune_sub_labels[immune_sub_labels=='gd T cells / NK-like cell (1)'] = "gd T cells\n/NK-like cell (1)"
immune_sub_labels[immune_sub_labels=="LSEC & Non-Inf Mac (9)"] = "LSEC &\nNon-Inf Mac (9)"
immune_sub_labels[immune_sub_labels=="Hep & non-Inf Mac (6)"] = "Hep &\nNon-Inf Mac (6)"
#### based on sample-type ##### 
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(set_1_cl_ord) ) 
counts$cluster= factor(counts$cluster, levels = as.character(set1_info$label[1:14]) ) 
counts$cluster= factor(counts$cluster, levels = as.character(set_2_cl_ord) ) 

counts$gcellType = gsub('\\(.*', '', as.character(counts$cluster))
counts$gcellType2 = counts$gcellType
counts$gcellType2[grepl(pattern = 'T cell', counts$gcellType )] = 'T cell'
counts$gcellType2[grepl(pattern = 'T cell', counts$gcellType )] = 'T cell &\nNK-like cell'
counts$gcellType2[grepl(pattern = 'Mac', counts$gcellType )] = 'Macrophage'
counts$gcellType2[grepl(pattern = 'B cell', counts$gcellType )] = 'B cell'
table(counts$gcellType2)

counts$strain=sapply(strsplit(counts$sample_type,'-') ,'[[', 1)
#write.csv(counts, 'figure_panel/set1-counts.csv')


ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  scale_fill_manual(values=mycolors)+
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

counts_norm$cluster2 = gsub(x =counts_norm$cluster, 'Inflammatory', 'Inf')
counts_norm$cluster2 = gsub(x =counts_norm$cluster2, 'and ', '& ')
ggplot(data=counts_norm, aes(x=cluster2, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Fraction of sample per cell type (%)')+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  #scale_fill_brewer(palette = "Set3")+
  #scale_fill_manual(values=mycolors)+
  #scale_x_discrete(labels = immune_sub_labels) #scales::wrap_format(9)
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=12.5,angle=90,color='black'),
        legend.title = element_blank()) +  
  xlab('')
  

 



########### dodge barplots -->>> why can't I add borders????
#counts2 = data.frame(read.csv('figure_panel/set2-counts.csv'))

gcelltype_df = aggregate(counts$Freq, by=list(gcellType=counts$gcellType2), FUN=sum)
gcelltype_ord = gcelltype_df$gcellType[order(gcelltype_df$x, decreasing = T)]
counts$gcellType2 = factor(counts$gcellType2, levels = as.character(gcelltype_ord) ) 
#aggregate(counts$Freq, by=list(gcellType=counts$test), FUN=sum)

ggplot(data=counts, aes(x=gcellType2, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity")+
  theme_classic()+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  #scale_fill_manual(values=mycolors)+
  theme(text = element_text(size=15),
        legend.title = element_blank(), axis.text.x = element_text(color='black',size=12, angle=90))+ 
  xlab('')+ylab('Counts')
  #scale_x_discrete(labels = c("Hep","LSEC", "Inf Mac", "Non-Inf\nMac","Stellate", "NK-like\n&T cells"))
  #scale_x_discrete(labels = c("B cell ","Cd3\nT cell","Erythroid","Hep","Inf Mac", "LSEC",
  #                            "Mac\n& DC", "Mature\nB cell","NK cell", "Non-Inf\nMac", "Stellate"))  
  scale_x_discrete(labels = c("Hep","LSEC","Erythroid", "NK cell",  "Mac\n& DC", "Cd3\nT cell",
                              "Inf Mac","Mature\nB cell","Non-Inf\nMac", "Stellate","B cell"))  

###### set-2 updated labels
set_2_cl_ord = c( "Hep (1)", "Hep (3)", "Hep (9)", "Hep (14)", "Hep (15)", "Hep (0)" , "Hep (2)",  "LSEC (6)" , 
                  "LSEC (4)", "Stellate (16)", "Non-Inflammatory Mac (8)",   "Non-Inflammatory Mac (13)", 
                  "Inflammatory Mac (11)",  "gd T cell (7)" , "Naive T cell (10)", "Mature B cell (12)", 
                  "pDC (17)", "Erythroid (5)")
#### immune subclustering
set_2_cl_ord = c("B cell & pDC (7)" , "B cell (3)", "NK-like cell (8)", "NK-like cell (1)" , "gd T cell (12)" , "Cd3 T cell (0)",
                 "Mac & LSEC (9)"  ,  "Mac & Hep (6)" ,  "Non-Inf Mac (2)" , "Non-Inf Mac (13)" , "Macrophage (11)"  ,
                 "Inf Mac (5)"  , "Inf Mac (4)", "Inf Mac (10)" ) 
      

###### set-1 updated labels
set_2_cl_ord = set_1_cl_ord                 
############## Dot plots ###############
############### Set-1 ########
merged_samples2 = merged_samples

set1_info_ord = set1_info[match(set_2_cl_ord, set1_info$label),]
matrix = as.matrix(set1_info_ord[,4:8])
markers = c()
for(i in 1:nrow(set1_info_ord)) markers = c(markers,matrix[i,])
markers <- markers[markers!='']
markers = unique(unname(markers))
div_num = 42#34 #58 #64
markers_final = unique(c(markers[1:div_num],'Ptprc', 'Cd68', markers[(div_num+1):length(markers)]))
markers_final = unique(c('Ptprc', 'Cd68',markers))
markers_final = unique(c('Ptprc', markers[1:div_num], 'Cd68', markers[(div_num+1):length(markers)]))

Idents(merged_samples2) <- factor(df_umap$label, levels = set_2_cl_ord)
DotPlot(merged_samples2, features = markers_final) + RotatedAxis()+ xlab('Markers')+
  ylab('')+theme(axis.text.x = element_text(size=12))
#split.by = "groups"


############## Varimax-related plots ##############
rot_data <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rot_data <- readRDS('Results/new_samples/varimax_rotated_object_new.rds') ## MT-removed
rot_data <- readRDS('Results/new_samples/immune_varimax_results.rds') ## set2-immune sub-population

rotatedLoadings <- rot_data$rotLoadings
Var_5_Load = data.frame(genes=rownames(rotatedLoadings),loading=round(rotatedLoadings[,5],3))
Var_15_Load = data.frame(genes=rownames(rotatedLoadings),loading=rotatedLoadings[,15])

head(Var_5_Load[order(Var_5_Load$loading, decreasing = T),],10)

scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
embedd_df_rotated <- data.frame(scores)

pc_num = 5
rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                     emb_val=embedd_df_rotated[,pc_num],
                     cluster=as.character(merged_samples$cluster),
                     Library_size=merged_samples$nCount_RNA,
                     num_expressed_genes=merged_samples$nFeature_RNA,
                     strain=merged_samples$strain,
                     label=factor(df_umap$label, levels = as.character(set_1_cl_ord)),
                     sample_name = merged_samples$sample_name
                     )

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=cluster))+geom_point(alpha=0.7,size=4)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(name='Immune\nsubcluster',values = mycolors)+theme_classic()+
  theme(text = element_text(size=22))+#, legend.title = element_blank()
  ggtitle(paste0('Set-2 Immune-subclusters over Varimax 1 and ', pc_num))


#### boxplots plots for varimax
ggplot(rot_df, aes(x=label, y=emb_val, fill=label))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = mycolors)+theme_classic()+
  +ylab(paste0('Varimax ', pc_num))+
  theme(text = element_text(size=18),
        axis.text.x = element_text(size=13,angle=90,color='black'),
        legend.title = element_blank(), legend.position = "none")+
  #scale_x_discrete(labels = xlabel2)+
  xlab('')

ggplot(rot_df, aes(x=label, y=emb_val, fill=strain))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  ylab(paste0('Varimax ', pc_num))+
  theme(text = element_text(size=18),
        axis.text.x = element_text(size=13,angle=90,color='black'),
        legend.title = element_blank(), legend.position = "none")+
  scale_x_discrete(labels = xlabel2)+xlab('')


df_umap$Varimax = abs(rot_df$emb_val)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=abs(Varimax)))+geom_point(alpha=0.5,size=1)+
  scale_color_viridis(direction = -1, option = 'inferno',name=paste0("Varimax ", pc_num))+theme_classic()+
  theme(text = element_text(size=22), legend.title=element_text(size=17))
#+ggtitle(paste0('Varimax ', pc_num, ' scores over set2\nimmune-subclusters UMAP'))
#'plasma', 'inferno', 'magma'

###### strain-specific varimax factors

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.5,size=1.5)+
  theme_classic()+ylab(paste0('varimax_',pc_num))+xlab("varimax 1")+
  scale_color_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  theme(text = element_text(size=16), legend.title = element_blank())# "#56B4E9"
  #+ggtitle(paste0('Distribution of cells based on strain\nover Varimax ', pc_num))
  
###### strain-specific varimax factors - sample

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point(alpha=0.5,size=1.4)+
  theme_classic()+ylab(paste0('Varimax-',pc_num))+xlab("Varimax-1")+theme_classic()+
  theme(text = element_text(size=18), legend.title = element_blank())# "#56B4E9"
#+ggtitle(paste0('Distribution of cells based on strain\nover Varimax ', pc_num))



data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(rot_df, aes(x=strain, y=emb_val, fill=strain))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  ylab(paste0('Varimax ', pc_num))+xlab('')+
  theme(text = element_text(size=17),axis.text.x = element_text(size=18,color='black'), legend.title = element_blank())+
  stat_summary(fun.data=data_summary)+stat_compare_means(label.x = 0.9, label.y = 3.8, size=5)

##########################################


###### Kmeans results on set-1 varimax-15 ###### 
kmeans_df <- readRDS('Results/kmeans_var15_1_oldsamples.rds')
ggplot(kmeans_df, aes(x=scores.Varimax_1,y=scores.Varimax_15,color=labels))+geom_point(alpha=0.7,size=4.5)+
  scale_color_manual(values=c("#E69F00","#56B4E9", "#999999"))+theme_classic()+
  theme(text = element_text(size=22), legend.title = element_blank())+# "#56B4E9"
  xlab('Varimax 1')+ylab('Varimax 15')
### check the lables on the total umap
umap_df2 <- data.frame(df_umap,
                       Kmeans_label=as.character(kmeans_df$labels))
ggplot(umap_df2, aes(UMAP_1, UMAP_2,color=Kmeans_label))+geom_point(alpha=0.4,size=3.5)+
  scale_color_manual(values=c("#E69F00","#56B4E9", "#999999"))+theme_classic()+
  theme(text = element_text(size=22), legend.title = element_blank())# "#56B4E9"



#####################################
################# Varimax correlation with technical factors
#####################################
merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

qc_variables <- c( "Library_size", "num_expressed_genes","mito_perc" )
main = 'Varimax-factors correlation with\ntechnical covariates (set-2)'
PCs_to_check <- c(1:25)
embedd_df_rotated_2 <- embedd_df_rotated[,PCs_to_check]
colnames(embedd_df_rotated_2) <- paste0('var-', PCs_to_check)
df.cov <- data.frame(embedd_df_rotated_2,
                     Library_size=merged_samples$nCount_RNA,
                     num_expressed_genes=merged_samples$nFeature_RNA, 
                     mito_perc=merged_samples$mito_perc)
df_cor <- cor(df.cov)
df_cor_sub <- df_cor[colnames(df_cor) %in% qc_variables, !colnames(df_cor)%in% qc_variables]
rownames(df_cor_sub) = c('Library size', '#Expressed\ngenes', 'Mt-gene\nfraction')
pheatmap(t(df_cor_sub), fontsize =10,fontsize_row=12,fontsize_col=12, main=main, 
         color = colorRampPalette(c("navy", "white", "red"))(50), )

pheatmap(t(df_cor_sub), fontsize =10,fontsize_row=12,fontsize_col=12, main=main, 
         color = colorRampPalette(c("blue3", "white", "red"))(50), )




##################################################
############# Random Forest models ################

result <- readRDS('Objects/cluster_pred_RFs_set1.rds')
result <- readRDS('Objects/cluster_pred_RFs_set2.rds')
RF_models <- result$models
preds <- result$preds
cm_list <- result$cm

cluster_num = 13
FeatImpDf <- data.frame(RF_models[[paste0('cluster_', cluster_num)]]$importance)
FeatImpDf$Varimax <-as.character(paste0('var ',1:50))
head(FeatImpDf)

FeatImpDf$Varimax <- factor(FeatImpDf$Varimax, levels=FeatImpDf$Varimax[order(FeatImpDf$MeanDecreaseGini)])
ggplot(FeatImpDf, aes(y=MeanDecreaseGini, x=as.factor(Varimax)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, fill='lightseagreen')+
  coord_flip()+theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 10),  
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15))+xlab('Varimax Factors')+ylab('Mean Decrease Gini Score')+
  ggtitle(paste0('Feature-importance score of a Random-Forest \n model tranied to predict cluster ',cluster_num))


FeatImpDf$Varimax <- factor(FeatImpDf$Varimax, 
                            levels=FeatImpDf$Varimax[order(FeatImpDf$MeanDecreaseGini, decreasing = T)])
FeatImpDf <- FeatImpDf[1:20,]
ggplot(FeatImpDf, aes(y=MeanDecreaseGini, x=as.factor(Varimax)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, fill='lightseagreen')+
  theme_classic()+#coord_flip()
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90),
        axis.text.y = element_text(color = "grey20", size = 10),  
        axis.title.x = element_text(color = "grey20", size = 13),
        axis.title.y = element_text(color = "grey20", size = 13))+
  xlab('Varimax Factors')+ylab('Mean Decrease Gini Score')+
  ggtitle(paste0('Top-20 feature-importance score of a Random-Forest\nmodel tranied to predict cluster ',cluster_num))



strain_pred_RF <- readRDS('Objects/strain_pred_RF_set1.rds')
strain_pred_RF <- readRDS('Objects/strain_pred_RF_set2.rds')


FeatImpDf <- data.frame(strain_pred_RF$importance)
FeatImpDf$Varimax <-as.character(paste0('Var ',1:50))
FeatImpDf$Varimax <- factor(FeatImpDf$Varimax, levels=FeatImpDf$Varimax[order(FeatImpDf$MeanDecreaseGini, decreasing = T)])
FeatImpDf <- FeatImpDf[1:20,]
head(FeatImpDf)
ggplot(FeatImpDf, aes(y=MeanDecreaseGini, x=as.factor(Varimax)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, fill='steelblue1')+
  theme_classic()+#coord_flip()
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90),
        axis.text.y = element_text(color = "grey20", size = 10),  
        axis.title.x = element_text(color = "grey20", size = 13),
        axis.title.y = element_text(color = "grey20", size = 13))+xlab('Varimax Factors')+ylab('Feature Importance\n(Mean Decrease Gini-Score)')
  #ggtitle('Top-20 feature-importance score of a Random-Forest \n model tranied to predict rat strains (set-1)')






############################################################
############### Quality control figures ###############
#############################################
# rat_DA_M09_WK_008_3pr_v3  rat_LEW_M09_WK_009_3pr_v3
# 'rat_DA_M_10WK_003' # 'rat_DA_01_reseq' 'rat_Lew_02'. rat_Lew_01 
INPUT_NAME = 'rat_LEW_M09_WK_009_3pr_v3'
MIT_CUT_OFF = 50
LIB_SIZE_CUT_OFF = 1000
NUM_GENES_DETECTED = 250
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
seur_raw[["mito_perc"]] <- PercentageFeatureSet(seur_raw, features = mito_genes_index)

df = data.frame(library_size= seur_raw$nCount_RNA, 
                mito_perc=seur_raw$mito_perc , 
                n_expressed=seur_raw$nFeature_RNA)


ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point(size=4,alpha=0.7)+
  scale_color_viridis(name='# expressed genes')+
  theme_classic()+xlab('Library size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red", size=0.9)+#labs(caption = INPUT_NAME)+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.9)+
  ggtitle('Library size and mitochondrial transcript cutoffs \nSet-2 Lew sample')+
  theme(text = element_text(size=18))


merged_samples <- readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
merged_samples <- readRDS('Objects/merged_samples_newSamples_MT-removed.rds')

df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2],
                      LibrarySize=merged_samples$nCount_RNA, 
                      numGenes=merged_samples$nFeature_RNA, 
                      MitoExp=merged_samples$mito_perc)

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=numGenes))+geom_point(alpha=0.7,size=3.5)+
  scale_color_viridis(name='Library size')+#Mt-gene\nfraction
  theme_classic()+theme(text = element_text(size=17))+
  ggtitle('Library size\nSet-2 map ')#Mitochondrial gene fraction


######################################
################ gProfiler results

gProfiler.df <- read.csv('figure_panel/gProfiler_csv/set2immune-var4-neg300.csv')
gProfiler.df <- read.csv('figure_panel/gProfiler_csv/set1-var15-Neg300.csv')
gProfiler.df <- read.csv('figure_panel/gProfiler_csv/set1-var15-pos300.csv')

gProfiler.df <- read.csv('figure_panel/gProfiler_csv/set1-var5-Neg300.csv')
gProfiler.df <- read.csv('figure_panel/gProfiler_csv/set1-var5-pos300.csv')

gProfiler.df <- gProfiler.df[,c('term_id','negative_log10_of_adjusted_p_value')]
colnames(gProfiler.df) = c('term', 'adj_pVal')
head(gProfiler.df)
dim(gProfiler.df)
summary(gProfiler.df$adj_pVal)
gProfiler.df <- gProfiler.df[1:15,]
gProfiler.df$term <- factor(gProfiler.df$term, levels=gProfiler.df$term[order(gProfiler.df$adj_pVal)])



##### visualization of the top TFs 
ggplot(gProfiler.df, aes(y=adj_pVal, x=term))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, aes(fill=adj_pVal))+# fill='lightseagreen'
  # scale_fill_gradientn(name='-log10(p.adj)',low='snow', high='blue')+
  scale_fill_gradientn(name='-log10(p.adj)',colours=brewer.pal(9,'Purples'),na.value = "transparent", #YlGn
                       breaks=c(2,4,6,8,10,12,15, 18, 20),
                       limits=c(1, 21),labels=c(2,'',6,'',10,'',15,'',20))+
  scale_x_discrete(limits = rev(levels(gProfiler.df$term[1:10])))+
  scale_y_continuous(limits = c(0, 21))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 15,angle = 90),
        axis.text.y = element_text(color = "grey20", size = 13),  
        axis.title.x = element_text(color = "grey20", size = 15),
        #legend.background = element_rect(fill="grey", size=0.5),
        axis.title.y = element_text(color = "grey20", size = 15))+
  ylab('-log10(adjusted p-value)')+
  xlab(paste0(''))#Transcription Factors

##### visualization of the top TFs - var-5
ggplot(gProfiler.df, aes(y=adj_pVal, x=term))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, aes(fill=adj_pVal))+# fill='lightseagreen'
  # scale_fill_gradientn(name='-log10(p.adj)',low='snow', high='blue')+
  scale_fill_gradientn(name='-log10(p.adj)',colours=brewer.pal(9,'Purples'),na.value = "transparent", #YlGn
                       breaks=c(2,4,6,8,10,12,15),
                       limits=c(1, 15),labels=c(2,'',6,'',10,'',15))+
  scale_x_discrete(limits = rev(levels(gProfiler.df$term[1:10])))+
  scale_y_continuous(limits = c(0, 15))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 15,angle = 90),
        axis.text.y = element_text(color = "grey20", size = 13),  
        axis.title.x = element_text(color = "grey20", size = 15),
        #legend.background = element_rect(fill="grey", size=0.5),
        axis.title.y = element_text(color = "grey20", size = 15))+
  ylab('-log10(adjusted p-value)')+
  xlab(paste0(''))#Transcription Factors


ggplot(gProfiler.df, aes(y=adj_pVal, x=as.factor(term)))+
  geom_bar(stat = "identity", colour = "black", width = 0.9, aes(fill=adj_pVal))+# fill='lightseagreen'
  scale_fill_gradientn(name='-log10(p.adj)',low='snow', high='blue')+
  coord_flip()+theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(color = "grey20", size = 13),  
        axis.title.x = element_text(color = "grey20", size = 15),
        #legend.background = element_rect(fill="grey", size=0.5),
        axis.title.y = element_text(color = "grey20", size = 15))+
  ylab('-log10(adjusted p-value)')+
  xlab(paste0('Transcription Factors'))


brewer.pal(9,'Blues')
display.brewer.pal(5,'Blues')


################################################### ############## 
############## implementaing a score to order the TF based on it ############## 
gProfiler.TF.DA <- read.csv('figure_panel/gProfiler_csv/set1-var15-Neg300.csv')
gProfiler.TF.lew <- read.csv('figure_panel/gProfiler_csv/set1-var15-pos300.csv')

gProfiler.TF.DA <- gProfiler.TF.DA[,c('term_id','negative_log10_of_adjusted_p_value')]
gProfiler.TF.lew <- gProfiler.TF.lew[,c('term_id','negative_log10_of_adjusted_p_value')]

colnames(gProfiler.TF.DA) = c('term', 'adj_pVal_DA')
colnames(gProfiler.TF.lew) = c('term', 'adj_pVal_LEW')

gProfiler.TF.lew$LEWscore = gProfiler.TF.lew$adj_pVal_LEW * (nrow(gProfiler.TF.lew):1)/nrow(gProfiler.TF.lew)
gProfiler.TF.DA$DAscore = gProfiler.TF.DA$adj_pVal_DA * (nrow(gProfiler.TF.DA):1)/nrow(gProfiler.TF.DA)

TF.merged = merge(gProfiler.TF.DA, gProfiler.TF.lew, by.x='term', by.y='term', all.x=T, all.y=T)
TF.merged2 = TF.merged[,c(1,3,5)]
TF.merged2[is.na(TF.merged2)] = 0
TF.merged2$dif = TF.merged2$LEWscore - TF.merged2$DAscore
TF.merged2.sort = TF.merged2[order(TF.merged2$dif,decreasing = T),]

TF.merged.vis = TF.merged2.sort[,2:3]
rownames(TF.merged.vis) = TF.merged2.sort$term
TF.merged.vis[TF.merged.vis==0]=NA
colnames(TF.merged.vis) = c('DA', 'LEW')
pheatmap(t(TF.merged.vis),cluster_rows = F, cluster_cols = F, color=inferno(40,direction = -1), fontsize_row  = 13, fontsize_col  = 10)
# set1-var15Pos-chea


#####################################
############## Ucell plots

library(stats)
library(ggpubr)
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))

kmeans_df <- readRDS('Results/kmeans_var15_1_oldsamples.rds')
high_score_cells = rownames(scores)[kmeans_df$labels %in% c('1','2')]


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
  ylab('POSITIVE REGULATION OF LYMPHOCYTE MIGRATION')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"

p2=ggplot(scores_2, aes(y=LYMPHOCYTE.MIGRATION.GOBP.GO.0072676_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('LYMPHOCYTE MIGRATION')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"

p3=ggplot(scores_2, aes(y=REGULATION.OF.LYMPHOCYTE.MIGRATION.GOBP.GO.2000401_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('REGULATION OF LYMPHOCYTE MIGRATION')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
gridExtra::grid.arrange(p1,p2,p3,nrow=1,ncol=3)



query_genesets = names(geneSets)[grepl(pattern = 'response to interferon-beta', names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features=geneIds(geneSets[names(geneSets) %in% query_genesets])))
colnames(scores_2) <- c('CELLULAR.RESPONSE.TO.INTERFERON.BETA', 'RESPONSE.TO.INTERFERON.BETA')
scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])

p1=ggplot(scores_2, aes(y=CELLULAR.RESPONSE.TO.INTERFERON.BETA, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('CELLULAR RESPONSE TO INTERFERON BETA')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
p2=ggplot(scores_2, aes(y=RESPONSE.TO.INTERFERON.BETA, x=strain))+geom_boxplot(aes(fill=strain))+
  ylab('RESPONSE TO INTERFERON BETA')+theme_classic()+
  scale_fill_manual(values=c("#E69F00", "#009E73"))+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
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




a_gene_set_name <- 'STAT4' # 'HNF4A' # GATA1
query_genesets = names(geneSets)[grepl(pattern = a_gene_set_name, names(geneSets),ignore.case = T)]

scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features= geneSets[names(geneSets) %in% query_genesets]))
colnames(scores_2) <- sapply(strsplit(colnames(scores_2),'\\.'),function(x) paste(x[1],x[5],x[7],sep = '_'))

scores_2$strain = sapply(strsplit(high_score_cells,'_'), function(x) x[2])
head(scores_2)

ggplot(scores_2, aes(y=STAT4_TH1_Mouse_UCell, x=strain))+geom_boxplot(aes(fill=strain))+scale_fill_brewer(palette = 'Dark2')+stat_compare_means()



########## evaluating the DA and LEW redundant gene sets enrichment in each

var_file = 'set1-var15-pos300'
gProfiler.df <- read.csv(paste0('figure_panel/gProfiler_csv/',var_file,'.csv'))

num_sets = nrow(gProfiler.df)
if(nrow(gProfiler.df)>25) num_sets = 15

candidateGenes.df.list = sapply(1:num_sets, function(i){
  candidateGenes = unlist(str_split(gProfiler.df$intersections[i], pattern = ','))
  candidateGenes.df = .getMapped_hs2model_df(ensembl, candidateGenes, model_animal_name)
  return(candidateGenes.df)
  },simplify = F
)
names(candidateGenes.df.list) = term_id_vec[num_sets]

pdf(paste0('figure_panel/gProfiler_csv/',var_file,'-uniqueGenes.pdf'))
for(i in 1:num_sets){
  print(i)
  term_id_vec = gProfiler.df$term_id[i]
  candidateGenes.df = candidateGenes.df.list[[i]]
  candidateGenes.rat = candidateGenes.df$rnorvegicus_homolog_associated_gene_name[candidateGenes.df$rnorvegicus_homolog_orthology_type == 'ortholog_one2one' & 
                                                                                    !is.na(candidateGenes.df$rnorvegicus_homolog_orthology_type)]
  
  # candidateGenes.rat = .getMapped_hs2model_df(ensembl,pos_specific, model_animal_name)
  scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix_sub, features= list(sig=c(unique(candidateGenes.rat)))))
  scores_2$strain = sapply(strsplit(high_score_cells,'_'), '[[', 2)
 
  p=ggplot(scores_2, aes(y=sig_UCell, x=strain))+geom_boxplot(aes(fill=strain))+
    scale_fill_brewer(palette = 'Dark2')+stat_compare_means()+ggtitle(paste0(term_id_vec, ' ', var_file))
  print(p)
}
dev.off()





#### check which gene-sets are enriched in both ways and identify the genes that are inconsitent between them

geneSets.var5.neg <- read.csv('figure_panel/gProfiler_csv/set1-var5-Neg300.csv')
geneSets.var5.pos <- read.csv('figure_panel/gProfiler_csv/set1-var5-pos300.csv')

geneSets.var15.neg <- read.csv('figure_panel/gProfiler_csv/set1-var15-Neg300.csv')
geneSets.var15.pos <- read.csv('figure_panel/gProfiler_csv/set1-var15-pos300.csv')

head(geneSets.var15.neg$term_id)
head(geneSets.var15.pos$term_id)

geneSets.var15.neg$term_id[1:15][!geneSets.var15.neg$term_id[1:15] %in% geneSets.var15.pos$term_id]
geneSets.var15.pos$term_id[1:15][!geneSets.var15.pos$term_id[1:15] %in% geneSets.var15.neg$term_id]

top_num = 20
same.geneset = geneSets.var15.neg$term_id[geneSets.var15.neg$term_id %in% geneSets.var15.pos$term_id]
same.geneset.top20 = geneSets.var15.neg$term_id[1:top_num][geneSets.var15.neg$term_id[1:top_num] %in% geneSets.var15.pos$term_id[1:top_num]]

# top 10: "SPI1"  "MECOM" "FLI1"  "TAL1"  "MYB"  
# top 20: "SPI1"  "RUNX1" "MECOM" "FLI1"  "TAL1"  "MYB"   "EGR1"  "RELA"  "CUX1" 
# top 30: "SPI1"  "RUNX1" "MECOM" "MITF"  "FLI1"  "TAL1"  "MYB"   "EGR1"  "RELA"  
###       "PPARG" "KLF1"  "GATA2" "CUX1"  "E2F1"  "CLOCK" "SOX2"  "STAT3" "STAT4" "CEBPB" "LMO2"  "DMRT1"

i = 1
same.geneset.top20[i]

genes.neg = geneSets.var15.neg$intersections[which(geneSets.var15.neg$term_id == same.geneset.top20[i])]
genes.pos = geneSets.var15.pos$intersections[which(geneSets.var15.pos$term_id == same.geneset.top20[i])]
genes.neg = unlist(str_split(genes.neg, ','))
genes.pos = unlist(str_split(genes.pos, ','))
neg_specific  = genes.neg[!genes.neg %in% genes.pos]
pos_specific  = genes.pos[!genes.pos %in% genes.neg]

var15.pos_specific.orth = .getMapped_hs2model_df(ensembl, pos_specific, model_animal_name)
var15.pos_specific.orth.genes = getUnemptyList(var15.pos_specific.orth$rnorvegicus_homolog_associated_gene_name)

var15.neg_specific.orth = .getMapped_hs2model_df(ensembl, neg_specific, model_animal_name)
var15.neg_specific.orth.genes = getUnemptyList(var15.neg_specific.orth$rnorvegicus_homolog_associated_gene_name)

pdf('Plots/var15.SPI1_pos_specific.orth.pdf')
for(i in 1:nrow(var15.pos_specific.orth)){
  gene2test = var15.pos_specific.orth.genes[i]
  if(!gene2test %in% rownames(my.matrix_sub)) next
  geneTest.df=data.frame(geneExp=my.matrix_sub[gene2test,],strain = sapply(strsplit(high_score_cells,'_'), '[[', 2))
  p=ggplot(geneTest.df, aes(x=strain,y=geneExp))+geom_boxplot()+ggtitle(gene2test)
  print(p)
}
dev.off()

pdf('Plots/var15.SPI1_neg_specific.orth.pdf')
for(i in 1:nrow(var15.neg_specific.orth)){
  gene2test = var15.neg_specific.orth.genes[i]
  if(!gene2test %in% rownames(my.matrix_sub)) next
  geneTest.df=data.frame(geneExp=my.matrix_sub[gene2test,],
                         strain = sapply(strsplit(high_score_cells,'_'), '[[', 2))
  p=ggplot(geneTest.df, aes(x=strain,y=geneExp))+geom_boxplot()+ggtitle(gene2test)
  print(p)
}
dev.off()


###### checking if the genes are redundant
var_num= 15 
varimax_res_set1 <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
var.load.df = data.frame(genes=rownames(varimax_res_set1$rotLoadings),
                         loading=varimax_res_set1$rotLoadings[,var_num])
var.load.df = var.load.df[order(var.load.df$loading, decreasing = T),]
head(var.load.df)
var15.pos300_genes = var.load.df$genes[1:300]
var15.neg300_genes = var.load.df$genes[nrow(var.load.df):(nrow(var.load.df)-300)]

sum(var15.pos300_genes %in% var15.neg300_genes)



###########################
###### Generating the total QC plot

#### new syncronized thresholds: 
MIT_CUT_OFF = 40
LIB_SIZE_CUT_OFF = 1500
NUM_GENES_DETECTED = 250

###### 
input_from_10x = paste0("Data/",c('rat_DA_M_10WK_003','rat_DA_01_reseq','rat_Lew_01', 'rat_Lew_02'),'/')
seur_raw <- lapply(input_from_10x, function(an_input_from_10x) CreateSeuratObject(counts=Read10X(an_input_from_10x, gene.column = 2),
                                          min.cells=0,min.features=1, 
                                          project = "snRNAseq"))

seur_raw.merge = merge(seur_raw[[1]], c(seur_raw[[2]], seur_raw[[3]], seur_raw[[4]]), # 
      add.cell.ids = c('rat_DA_M_10WK_003','rat_DA_01_reseq','rat_Lew_01', 'rat_Lew_02'), 
      project = "rat_data", 
      merge.data = TRUE)
dim(seur_raw.merge)
rm(list=ls())

seur_genes_df <- read.delim(paste0(input_from_10x[1],'features.tsv.gz'), header = F)
seur_raw.merge[['RNA']] <- AddMetaData(seur_raw.merge[['RNA']], seur_genes_df$V2, col.name = 'symbol')
seur_raw.merge[['RNA']] <- AddMetaData(seur_raw.merge[['RNA']], seur_genes_df$V1, col.name = 'ensembl')

MIT_PATTERN = '^Mt-'
mito_genes_index <- grep(pattern = MIT_PATTERN, seur_raw.merge[['RNA']]@meta.features$symbol )
seur_raw.merge[['RNA']]@meta.features$symbol[mito_genes_index]
seur_raw.merge[["mito_perc"]] <- PercentageFeatureSet(seur_raw.merge, features = mito_genes_index)
summary(seur_raw.merge[["mito_perc"]]$mito_perc )
libSize <- colSums(GetAssayData(seur_raw.merge@assays$RNA))

seur_raw = seur_raw.merge 
df = data.frame(library_size= seur_raw$nCount_RNA, 
                mito_perc=seur_raw$mito_perc , 
                n_expressed=seur_raw$nFeature_RNA)


## Visualization of QC metrics
ggplot(df, aes(x=library_size, y=mito_perc, color=n_expressed))+geom_point(alpha=0.7)+
  scale_color_viridis(name='#expressed\ngenes')+
  theme_classic()+xlab('Library Size')+ylab('Mitochondrial transcript percent')+
  geom_hline(yintercept= MIT_CUT_OFF, linetype="dashed", color = "red")+
  geom_vline(xintercept = LIB_SIZE_CUT_OFF, linetype="dashed", color = "red3", size=0.5)+
  theme(text = element_text(size=16))


#######################################################
##### saving varimax loadings as supp material  #######
#######################################################

rot_data <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rot_data <- readRDS('~/RatLiver/Results/new_samples/varimax_rotated_object_new.rds') ## MT-removed


loadings_list = sapply(1:ncol(rot_data$rotLoadings), function(i)rot_data$rotLoadings[,i], simplify = F)
loadings_df = data.frame(do.call(cbind, loadings_list))
colnames(loadings_df) = paste0('varPC_', 1:ncol(loadings_df))
head(loadings_df)

#write.csv(loadings_df, '~/RatLiver/Results/old_samples/varimax_loadings_TLH.csv')
write.csv(loadings_df, '~/RatLiver/Results/new_samples/varimax_loadings_ImmuneEnriched.csv')





