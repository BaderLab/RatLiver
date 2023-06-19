source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)
########################################################
############ Importing the nuc-seq data  #########

##### snnRNA-seq data
merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
# merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')

####### scRNA-seq TLH and immune-enriched maps
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

load(new_data_scCLustViz_object)
load(old_data_scClustViz_object)

####### scRNA-seq TLH and snRNA-seq map integration
sham_sn_merged_scCLustViz_object =  '~/rat_sham_sn_data/standardQC_results/TLH_snMap_merged.RData' ### find the results on run1
load(sham_sn_merged_scCLustViz_object)


merged_samples = your_scRNAseq_data_object
merged_samples = merged_data
##########################################################################
#################### multiple gene expression over UMAP #################### 
##########################################################################

gene_list = c("Alb", "Apoa1", "Apoc1", "Apoc3", "Apoe", "Fabp1", "Itih4", "Orm1", "Pigr",
              "Serpina1", "Tf", "Ttr", "Sds", "Pck1", "Arg1", 
              "Ahr", "Akr1c1", "Cyp27a1", "Cyp7a1", "Glul", "Notum", "Rcan1", 
              "Cyp2e1", "Hamp") #"Cyp2fa1", "Mt2a"

gene_list = c("Cyp1a2", "Aldh6a1", "Cyp8b1", "Cyp7a1", "Fasn", ## set1
              "Acly", "Agt", "Aldob", "Arg1", "Cps1", "Cyp2a1", "Insig1", "Pck1", "Sds", "Ttr2", "Arg1", ## set2
              "Fth1", "G6pc", "Slco1b2", "Igfbp1", "Insig1", "Scd", "Itih3", "Fgb", "C3", "Calr", "Dbi", "Apoa2", "Cox8a", "Rpl10", ## set3
              "Ctsl", "Fcgr2b", "Sparc", "Stab2", "Ramp2", "Calcrl", "Ifi27", "Fam167b", "Eng", "Id3", "Aqp1", "Gpr182", "Lyve1", "Clec14a", ## set4
              "Aqp1", "Epcam", "Sox9", "Anpep", "Anxa4", ## set5
              "Marco", "Vsig4", "Clec4f", "Cd5l", "Hmox1", "Aif1", "Cd68", "Ptprc", "Cd163", "C1qb", "C1qc", "C1qa", "Slc11a", "Ptas1", "Ccl6", ## set6
              "Klrd1", "Gzmk", "Cd7", "Ccl5", "Gzma", "Nkg7", "Tmsb4x", "Il2rb", ## set7
              "Sparc", "Calcrl", "Ecm1", "Igfbp7", "Rbp1", "Col3a1", "Lyve1", "Bgn", "Colec10", "Colec11", ## set8
              "Rspo3", "Clec4g", "Wnt2", "Wnt9b", "Thad", "Thbd", "Fcn2") ## set9

gene_list = c("Cd74", "Lyz2", "Rt1-db1", "Rt1-da", "S100a11", "Ptprc", "RT1-Ba", "RT1-Bb", "RT1-Da", "RT1-Db1")

length(gene_list)
gene_list = gene_list[gene_list %in% rownames(merged_samples)]
length(gene_list)

merged_samples
embedding_name = 'umap_h'
embedding_name = 'umap'

dev.off()
pdf('snRNAseq_marker_umaps.pdf')
#pdf('ImmuneEnriched_marker_umaps.pdf')
pdf('TLH_marker_umaps.pdf')
pdf('sc_snRNAseq_integrated_marker_umaps.pdf')

for(i in 1:length(gene_list)){
  gene_name = gene_list[i]  #Sirpa
  rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset
  
  df_umap <- data.frame(UMAP_1=getEmb(merged_samples, embedding_name)[,1], 
                        UMAP_2=getEmb(merged_samples, embedding_name)[,2], 
                        expression_val=GetAssayData(merged_samples)[gene_name,])
  
  
  p=ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
    geom_point(alpha=0.6,size=1.2)+
    scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
    theme(text = element_text(size=16.5),
          plot.title = element_text(hjust = 0.5),
          legend.title=element_text(size = 10))+
    ggtitle(paste0(gene_name)) 
  print(p)
}
dev.off()
rm(merged_samples)


##########################################################################
#################### Single gene expression over UMAP #################### 
##########################################################################
Resolution = 2.5
merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)
table(merged_samples$SCT_snn_res.2.5)
merged_samples$clusters = merged_samples$SCT_snn_res.2.5


gene_name = 'Cd5l'  #Sirpa
rownames(merged_samples)[grep(gene_name, rownames(merged_samples))] ## check if the gene is present in the dataset
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap_h')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap_h')[,2], 
                      clusters=merged_samples$clusters,
                      expression_val=GetAssayData(merged_samples)[gene_name,],
                      umi=colnames(merged_samples))



ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=expression_val))+
  geom_point(alpha=0.6,size=1.2)+
  scale_color_viridis('Expression\nValue', direction = -1)+theme_classic()+
  theme(text = element_text(size=16.5),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_text(size = 10))+
  ggtitle(paste0(gene_name)) 



