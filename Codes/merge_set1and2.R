source('Codes/Functions.R')
Initialize()
library(RColorBrewer)


############ set-1 map
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
cluster_df = data.frame(umi=colnames(your_scRNAseq_data_object),
                        clusters=as.character(sCVdata_list$res.0.6@Clusters))
set1_info <- read.csv('figure_panel/set-1-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info$clusters = as.character(set1_info$clusters)
cluster_df_set1 = merge(cluster_df_set1, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
cluster_df_set1 <- cluster_df_set1[match(colnames(your_scRNAseq_data_object),cluster_df_set1$umi),]
cluster_df_set1$umi == colnames(your_scRNAseq_data_object)
cluster_df_set1$sample_name = unlist(lapply(str_split(cluster_df_set1$umi, pattern = '_'), 
                            function(x) paste0(x[-length(x)],collapse = '_')))
head(cluster_df_set1)
table(cluster_df_set1$sample_name)

############ set-2 map
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
load(new_data_scCLustViz_object)
cluster_df_set2 = data.frame(umi=colnames(your_scRNAseq_data_object),
                        clusters=as.character(sCVdata_list$res.0.6@Clusters))
set2_info <- read.csv('figure_panel/set-2-final-info.csv')
set2_info <- set2_info[1:17,]
colnames(set2_info)[1] = 'clusters'
set2_info$clusters = as.character(set2_info$clusters)
cluster_df_set2 = merge(cluster_df_set2, set2_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
cluster_df_set2 <- cluster_df_set2[match(colnames(your_scRNAseq_data_object),cluster_df_set2$umi),]
cluster_df_set2$umi == colnames(your_scRNAseq_data_object)
cluster_df_set2$sample_name = unlist(lapply(str_split(cluster_df_set2$umi, pattern = '_'), 
                                            function(x) paste0(x[-length(x)],collapse = '_')))
head(cluster_df_set2)


############ merging the umi/cluster info of set1 and set2 
cluster_label_df = rbind(cluster_df_set1, cluster_df_set2)
head(cluster_label_df)
cluster_label_df[(!cluster_label_df$umi %in% umap_df$umi),]
table(cluster_label_df$sample_name)
cluster_label_df$sample_name = ifelse(cluster_label_df$sample_name=='rat_DA_01_reseq', 'rat_DA_01',cluster_label_df$sample_name)


set1and2merged <- readRDS('~/XSpecies/Results/preproc_rats/merged/all_merged_samples_seur.rds')
table(set1and2merged$sample_name)

umap_df = data.frame(map=ifelse(set1and2merged$data_label=='old', 'set1', 'set2'), 
                     cluster=set1and2merged$clusters,
                     strain=set1and2merged$strain, 
                     Embeddings(set1and2merged, reduction = 'umap'),
                     umi=gsub('merged_','',colnames(set1and2merged)))
umap_df$umi = paste0(umap_df$umi, '-1')
### adding the cell labels and cluster into to the merged dataset
umap_df_label = merge(umap_df, cluster_label_df, by.x='umi', by.y='umi', all.x=T, all.y=F) 
umap_df_label$label2 = gsub("\\s*\\([^\\)]+\\)","",as.character(umap_df_label$label)) ### making the cell type annotations more general
umap_df_label$cluster = gsub('cluster_', '', umap_df_label$cluster) ### general clustering of the merged maps
umap_df_label$sample = unlist(lapply(str_split(umap_df_label$umi, pattern = '_'), 
                                            function(x) paste0(x[-length(x)],collapse = '_')))

sum(is.na(umap_df_label$label)) ### 5190 cells have been removed in the newer analysis (QC updates?)
table(umap_df_label$map[is.na(umap_df_label$label)])

# Define the number of colors you want
nb.cols <- length(names(table(umap_df_label$cluster))) #18 #14
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)


##### making UMAP plots ####
ggplot(umap_df_label, aes(x=UMAP_1, y=UMAP_2, color=label2))+geom_point(size=3,alpha=0.6)+
  theme_classic()+ scale_color_manual(values = mycolors)
ggplot(umap_df_label, aes(x=UMAP_1, y=UMAP_2, color=map))+geom_point(size=1.32,alpha=0.3)+
  theme_classic()
ggplot(umap_df_label, aes(x=UMAP_1, y=UMAP_2, color=cluster))+geom_point(size=3,alpha=0.9)+
  theme_classic()+ scale_color_manual(values = mycolors)


#### barplots of the contribution of each map to general clusters  ##### 
counts <- ddply(umap_df_label, .(umap_df_label$map, umap_df_label$cluster), nrow)
names(counts) <- c("map", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(0:20) ) 


ggplot(data=counts, aes(x=cluster, y=Freq, fill=map)) +
  geom_bar(stat="identity",color='black')+theme_classic()+
  ylab('Counts')+xlab('Clusters')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+xlab('')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=map)) +
  geom_bar(stat="identity",color='black')+theme_classic()+
  ylab('Contribution of map per cluster (%)')+xlab('Cluster')+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=11,angle=90,color='black'),
        legend.title = element_blank()) +  xlab('') 



###### evaluating the varimax results ######
varimax_res = readRDS('~/XSpecies/Results/preproc_rats/merged/varimax_results_All_merged_rats.rds')
embedd_df_rotated = data.frame(varimax_res$rotScores)
colnames(embedd_df_rotated) = paste0('Varimax_', 1:ncol(embedd_df_rotated))
PCs_to_check = ncol(embedd_df_rotated)

pdf('Plots/VarimaxPCA_set1set2_merged.pdf',width = 14,height=21) 

for(i in 1:PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       clusters=umap_df$cluster,
                       sample = umap_df_label$sample,
                       strain = umap_df$strain,
                       label=umap_df_label$label2)
  
  p1=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=label))+geom_point(alpha=0.8, size=1.1)+
    theme_classic()+ylab(paste0('Varimax ',i))+scale_color_manual(values=mycolors)+
    theme(text = element_text(size=15),legend.title = element_blank())#
  
  p2=ggplot(rot_df, aes(x=label, y=emb_val, fill=label))+geom_boxplot()+
    theme_classic()+scale_fill_manual(values=mycolors)+
    theme(text = element_text(size=18),
          axis.text.x = element_text(size=11,angle=90,color='black'),
          legend.title = element_blank(), legend.position = "none")+xlab('')+ylab(paste0('Varimax ',i))
  
  p3=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample))+geom_point(alpha=0.6,size=1.5)+
    theme_classic()+ylab(paste0('Varimax ',i))+
    #scale_color_manual(values = c("pink1", "skyblue1"))+
    theme(text = element_text(size=15),legend.title = element_blank())
  
  p4=ggplot(rot_df, aes(x=sample, y=emb_val, fill=sample))+geom_boxplot()+
    theme_classic()+
    #scale_fill_manual(values = c("pink1", "skyblue1"))+
    theme(text = element_text(size=10,angle = 45),legend.position='none')+
    ylab(paste0('Varimax ',i))+xlab('Strain')
  
  p5=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.6,size=1.5)+
    theme_classic()+ylab(paste0('Varimax ',i))+
    scale_color_manual(values = c("pink1", "skyblue1"))+
    theme(text = element_text(size=15),legend.title = element_blank())
  
  
  p6=ggplot(rot_df, aes(x=strain, y=emb_val, fill=strain))+geom_boxplot()+
    theme_classic()+scale_fill_manual(values = c("pink1", "skyblue1"))+
    theme(text = element_text(size=15),legend.title = element_blank())+
    ylab(paste0('Varimax ',i))+xlab('Strain')
  
  umap_df_label$Varimax = rot_df$emb_val
  p7=ggplot(umap_df_label, aes(x=UMAP_1, y=UMAP_2, color=abs(Varimax)))+geom_point(alpha=0.5,size=1.1)+
    scale_color_viridis(direction = -1, option = 'inferno',name=paste0("Varimax-", pc_num))+theme_classic()+
    theme(text = element_text(size=22), legend.title=element_text(size=17))
  
  nb.cols <- length(names(table(umap_df_label$label2))) #18 #14
  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
  
  p8=ggplot(umap_df_label, aes(x=UMAP_1, y=UMAP_2, color=label2))+geom_point(alpha=0.8, size=1.2)+
    theme_classic()+scale_color_manual(name='clusters',values = mycolors)+
    theme(text = element_text(size=15),legend.title = element_blank())#
  
  gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=2,nrow=4)
}
dev.off()            



