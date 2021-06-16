source('Codes/Functions.R')
Initialize()

#### Adding immune-sub annotations
annot_df <- read.csv('cluster_labels/newSamples.csv') #newSamples.csv
annot_df$cluster <- as.character(0:(nrow(annot_df)-1))
head(annot_df)

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
new_data_scCLustViz_object_endothelial <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_EndothelialSub.RData"

load(new_data_scCLustViz_object)


umap_df=data.frame(Embeddings(your_scRNAseq_data_object, 'umap'),
                   cluster=as.character(sCVdata_list$res.0.6@Clusters))

ggplot(umap_df, aes(UMAP_1, UMAP_2, color=cluster))+geom_point()+theme_classic()+
  ggtitle('Set-2')+
  scale_color_manual(values = colorPalatte)+
  theme(legend.title=element_blank())

## old samples
annot_df$cluster2=c(rep('Hep',3),'LSEC','Hep', 'Mac', 'Hep', 'Stellate', 'Hep', 'Plasma', 
                    'Mac', 'LSEC', 'Hep', 'T cell', 'Stellate', rep('Hep',2))

saveRDS(annot_df, 'cluster_labels/old_labels.rds')

## new samples
annot_df$cluster2=c(rep('Hepatocytes',3), rep('Central Venous LSECs',2),'Erythroid Cells', 
  'Central Venous LSECs', 'Alpha-Beta/Gamma-Delta T Cells 1 and NK-like Cells',
  'Mac', 'Hepatocytes', 'Erythroid Cells', 'Inf Mac', 'Plasma', 'Mac', rep('Hepatocytes',2),
  'Hepatic Stellate Cells', 'B Cells, Non-Inf Mac')

annot_df$labels=c(rep('Hep',3), rep('LSEC',2),'Erythroid', 
                    'LSEC', 'T cell',
                    'Mac', 'Hep', 'Erythroid', 'Inf Mac', 'Plasma', 'Mac', rep('Hep',2),
                    'Stellate', 'Non-Inf Mac')

saveRDS(annot_df, 'cluster_labels/new_labels.rds')
head(umap_df)
umap_df_2 <- merge(umap_df, annot_df, 'cluster', 'cluster', all.x=T)
head(umap_df_2)
ggplot(umap_df_2, aes(UMAP_1, UMAP_2, color=cluster2))+geom_point()+theme_classic()+
  ggtitle('Set-2')+
  scale_color_manual(values = colorPalatte)+
  theme(legend.title=element_blank())
