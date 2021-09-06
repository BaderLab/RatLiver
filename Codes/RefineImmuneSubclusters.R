source('Codes/Functions.R')
source('Codes/FactorAnalysisUtils.R')
Initialize()


##### Functions #####
clean_name <- function(x){
  x=x[-length(x)] 
  if(x[1]=='merged') x=x[-1]
  return(paste(x, collapse = '_'))
}


qc_variables <- c( "Library_size", "num_expressed_genes","mito_perc" )

check_qc_cor <- function(merged_samples, pca_embedding_df, main){
  df_cor <- data.frame(pca_embedding_df,
                       Library_size=merged_samples$nCount_RNA,
                       num_expressed_genes=merged_samples$nFeature_RNA, 
                       mito_perc=merged_samples$mito_perc)
  df_cor <- cor(df_cor[,!colnames(df_cor)%in% c('clusters')])
  df_cor_sub <- df_cor[colnames(df_cor)%in%qc_variables, !colnames(df_cor)%in%qc_variables]
  pheatmap(df_cor_sub, fontsize_row=10, main=main, color = magma(20,direction = 1))
}



new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
load(new_data_scCLustViz_object_Immune)
old_immune_data = your_scRNAseq_data_object
old_immune_data$subclusters = paste0('old-', as.character(sCVdata_list$RNA_snn_res.1@Clusters))
  
new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData"
load(new_data_scCLustViz_object_Immune)
new_immune_data = seur
new_immune_data$subclusters = as.character(new_immune_data$SCT_snn_res.1)

df_old <- data.frame(umi=colnames(old_immune_data), old_subclust=old_immune_data$subclusters) 
df_new <- data.frame(umi=colnames(new_immune_data), new_subclust=new_immune_data$subclusters, Cluster=new_immune_data$final_cluster)
merged_df = merge(df_new, df_old, by.x='umi', by.y='umi', all.x=T)
pheatmap(table(paste0('new-',merged_df$new_subclust), merged_df$old_subclust))
pheatmap(table(paste0('new-',merged_df$new_subclust), merged_df$Cluster))




load('Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData')
merged_samples = seur
merged_samples$subclusters = sCVdata_list$res.1@Clusters

merged_samples = readRDS('Results/new_samples/Immune_subclusters_labelCorrect_HepRemoved.rds')
merged_samples = ScaleData(merged_samples)
merged_samples <- RunPCA(merged_samples, features = rownames(merged_samples))  
merged_samples$sample_name = merged_samples$strain

#### check if removal of mt-genes would help with the annotation
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      #clusters=paste0('', sCVdata_list$res.0.6@Clusters),
                      clusters=merged_samples$cluster,
                      #large_clusters=merged_samples$final_cluster,
                      sample=merged_samples$sample_name,
                      strain=merged_samples$strain,
                      Ptprc=GetAssayData(merged_samples)['Ptprc',], 
                      umi=colnames(merged_samples))

##### adding the final annotations to the dataframe
set1_info <- read.csv('figure_panel/set-2-immunesub-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info <- set1_info[1:14,]
set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)
df_umap$label = merged_samples$label

rot_data <- readRDS('Results/new_samples/immune_varimax_results_cl17Inc.rds')

loading_matrix = Loadings(merged_samples, 'pca')
gene_exp_matrix = GetAssayData(merged_samples, assay = 'RNA') #SCT
rot_data <- get_varimax_rotated(gene_exp_matrix, loading_matrix)

rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))


### varimax PCA
rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'Varimax on last 2 rat samples - immune cells')
ndims = ncol(rot_data$rotScores)
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)))
perc_variance_threshold = 0.5
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))
PCs_to_check <-  c(1:20,31:50)


embedd_df_rotated <- data.frame(scores)
embedd_df_rotated <- embedd_df_rotated[,PCs_to_check]
colnames(embedd_df_rotated) <- paste0('Varimax_', PCs_to_check)
check_qc_cor(merged_samples, embedd_df_rotated, 
             main='Varimax-PC embeddings correlation with technical covariates')

pheatmap(cor(embedd_df_rotated[,c(PCs_to_check)]))

xlabel = c("B cell (3)", "Cd3 T cell (0)", "DC & B cell (7)", "gd T cell (12)", "gd T cells\n/NK-like cell (1)",
           "Hep &\nNon-Inf Mac (6)", "Inf Mac (10)", "Inf Mac (4)" ,"Inf Mac (5)", "LSEC &\nNon-Inf Mac (9)", "NK-like cell (8)", 
           "Non-Inf Mac (11)", "Non-Inf Mac (13)", "Non-Inf Mac (2)") 


nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols) #Pastel1


##### Visualizing the PCA plots after the varimax rotation #####
plot_dir = 'Plots/'
pdf(paste0(plot_dir, 'VarimaxPCA_scTrans_mergedNewSamples_immuneCells_c17Inc_HepRemoved.pdf'),width = 14,height=17) 

for(i in PCs_to_check){ 
  pc_num = i
  rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       clusters=merged_samples$cluster,
                       #large_clusters=merged_samples$final_cluster,
                       Library_size=merged_samples$nCount_RNA,
                       num_expressed_genes=merged_samples$nFeature_RNA,
                       sample_name = merged_samples$strain,#sapply(str_split(merged_samples$sample_name, '_'), '[[', 2),
                       label=merged_samples$label)
  
  p1=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=label))+geom_point(alpha=0.8, size=1.1)+
    theme_classic()+ylab(paste0('Varimax ',i))+scale_color_manual(values=mycolors)+
    theme(text = element_text(size=15),legend.title = element_blank())#
  
  p2=ggplot(rot_df, aes(x=label, y=emb_val, fill=label))+geom_boxplot()+
    theme_classic()+scale_fill_manual(values=mycolors)+
    theme(text = element_text(size=18),
          axis.text.x = element_text(size=11,angle=90,color='black'),
          legend.title = element_blank(), legend.position = "none")+
    scale_x_discrete(labels = xlabel)+xlab('')+ylab(paste0('Varimax ',i))
  
  p3=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point(size=1.5)+
    theme_classic()+ylab(paste0('Varimax ',i))+
    scale_color_manual(values = c("pink1", "skyblue1"))+
    theme(text = element_text(size=15),legend.title = element_blank())
  
  p4=ggplot(rot_df, aes(x=sample_name, y=emb_val, fill=sample_name))+geom_boxplot()+
    theme_classic()+scale_fill_manual(values = c("pink1", "skyblue1"))+
    theme(text = element_text(size=15),legend.title = element_blank())+
    ylab(paste0('Varimax ',i))+xlab('Strain')
  
  df_umap$Varimax = rot_df$emb_val
  p5=ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=abs(Varimax)))+geom_point(alpha=0.5,size=1.1)+
    scale_color_viridis(direction = -1, option = 'inferno',name=paste0("Varimax-", pc_num))+theme_classic()+
    theme(text = element_text(size=22), legend.title=element_text(size=17))
  
  p6=ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label))+geom_point(alpha=0.8, size=1.2)+
    theme_classic()+scale_color_manual(name='clusters',values = mycolors)+
    theme(text = element_text(size=15),legend.title = element_blank())#
  
  gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2,nrow=3)
}
dev.off()                      

rot_df$sample_name[rot_df$label=="Hep & non-Inf Mac (6)"]

ggplot(rot_df, aes(x=label, y=emb_val, fill=sample_name))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=18),
        axis.text.x = element_text(size=11,angle=90,color='black'),
        legend.title = element_blank())+ #legend.position = "none"
  scale_x_discrete(labels = xlabel)+xlab('')+ylab(paste0('Varimax ',i))


var_PC = 6
loading_emd_df=data.frame(gene=rownames(loading_matrix),
                          loading=loading_matrix[,var_PC])

loading_emd_df <- loading_emd_df[order(loading_emd_df$loading, decreasing = T),]
rownames(loading_emd_df) = NULL
dev.off()
gridExtra::grid.table(head(loading_emd_df,20))
dev.off()
gridExtra::grid.table(tail(loading_emd_df,20))




############## Mac specific dop plot #############
## subcluster each of the clusters mentioned in the slides
## based on the initial assumptions on the subcluster identity >> check the expression of a set of above markers using dotplots
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 14 #18 #14
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included.RData"
load(new_data_scCLustViz_object_Immune)
merged_samples <- seur
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      clusters=as.character(merged_samples$SCT_snn_res.1),
                      Ptprc=GetAssayData(merged_samples)['Ptprc',],
                      Cd68=GetAssayData(merged_samples)['Cd68',],
                      umi=colnames(merged_samples))


set1_info <- read.csv('figure_panel/set-2-immunesub-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info <- set1_info[1:14,]
set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)


Mac_subclusters <- c(2, 4, 5, 6, 9, 10 , 11, 13)
df_umap$Mac <- df_umap$cluster %in% Mac_subclusters
df_umap$label_Mac <- ifelse(df_umap$Mac, df_umap$label, 'Other')
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=label_Mac))+geom_point(alpha=0.8, size=2.5)+
  theme_classic()+scale_color_manual(name='clusters',values = c(mycolors[c(6:10,12:14)],'grey86'))+
  theme(text = element_text(size=15),legend.title = element_blank())#


ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Cd68))+geom_point(alpha=0.7,size=1)+
  scale_color_viridis(direction = -1)+theme_classic()+
  theme(text = element_text(size=20))

set_2_cl_ord = c("Non-Inf Mac (13)", "Non-Inf Mac (2)", "Hep & non-Inf Mac (6)",
                 "LSEC & Non-Inf Mac (9)", "Non-Inf Mac (11)","Inf Mac (5)" ,"Inf Mac (10)", "Inf Mac (4)")
mac_markers <- read.csv('~/RatLiver/Results/new_samples/Mac_markers.csv')
table(mac_markers$Annotation)
mac_markers$Marker[!mac_markers$Marker %in% rownames(merged_samples)]

merged_samples_sub <- merged_samples[,df_umap$Mac]
merged_samples_sub$orig.ident <- df_umap$label_Mac[df_umap$Mac]
Idents(merged_samples_sub) <- factor(merged_samples_sub$orig.ident, levels = set_2_cl_ord)
DotPlot(merged_samples_sub, features = c('Cd68',mac_markers$Marker[mac_markers$Marker %in% rownames(merged_samples)])) + 
  RotatedAxis() #+ ggtitle(names(Cluster_markers_sorted)[i])



