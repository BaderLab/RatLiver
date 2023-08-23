source('Codes/Functions.R')
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')
library(Hmisc)
library(RColorBrewer)
library(viridis)
library(reshape2)
########################################################
############ Importing the nuc-seq data  #########

## importing the gene expression data
merged_samples_sub = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
merged_samples_all = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_allfeatures.rds')


'Cd68'%in% rownames(GetAssayData(merged_samples_all, 'scale.data'))
'Cd68'%in% rownames(GetAssayData(merged_samples_sub, 'scale.data'))

markers_test <- c('Alb', 'Apoa1', 'Apoc1', 'Apoc3', 'Apoe', 'Fabp1', 
                  'Itih4', 'Orm1', 'Pigr', 'Serpina1', 'Tf', 'Ttr')

markers_test %in% rownames(GetAssayData(merged_samples_all, 'scale.data'))
markers_test[!markers_test %in% rownames(GetAssayData(merged_samples_all, 'scale.data'))]

dim(merged_samples_all)
dim(merged_samples_sub)

#### checking if the order of cells have been preserved in the data containing all the genes 
sum(colnames(merged_samples_all) != colnames(merged_samples_sub))
head(merged_samples_all)
head(merged_samples_sub)

Resolution = 2.5
resolutions = Resolution
merged_samples_sub <- FindClusters(merged_samples_sub, resolution = Resolution, verbose = FALSE)
table(merged_samples_sub$SCT_snn_res.2.5)

merged_samples_all$SCT_snn_res.2.5 = merged_samples_sub$SCT_snn_res.2.5
merged_samples = merged_samples_all


####################################################
#Resolution = 2.5
#merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)

#Resolution = 0.6
#merged_samples <- FindClusters(merged_samples, resolution = Resolution, verbose = FALSE)

#mapping_table = table(merged_samples$SCT_snn_res.2.5, merged_samples$SCT_snn_res.0.6)

#### defining the mega clusters based on 0.6 resolution
Hep0 = as.character(c(4, 7, 10, 16, 25, 27, 0)) # 27 is the little tail
Hep1 = as.character(c(31, 6, 12, 1, 17, 15)) # 31 is the tail
Hep2 = as.character(c(21, 9, 2, 5)) #
Hep3 = as.character(c(23, 3, 8, 13, 32, 26)) 

##### updated version (V1) - July 14th
Hep0 = c(4, 7, 10, 16, 25, 27, 0)
Hep1 = c(31, 6, 12, 1, 17, 15, 26)
Hep2 = c(21, 9, 2, 5)
Hep3 = c(23, 3, 8, 13, 32)

########## updated version (V2) - July 19th
Hep3 = c("23", "3", "8", "13", "32")
Hep2 = c("21", "9", "2", "5")
Hep0 = c("4", "7", "10", "16", "25", "27", "0")
Hep1 = c("31", "6", "12", "1", "17", "15", "26")

##########  updated version (V3) - July 19th - FINAL PUBLICATION FIGURE
Hep3 = c("3", "8", "13", "23", "32")
Hep2 = c("21", "9", "2", "5")
Hep0 = c("25", "27", "4", "16", "7", "0", "10")
Hep1 = c("1", "6", "15", "17", "12", "31", "26")



merged_samples$clusters = as.character(merged_samples$SCT_snn_res.2.5)
merged_samples$Hep_clusters = ifelse(merged_samples$clusters %in% as.character(Hep0), 'Hep0', merged_samples$clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep1), 'Hep1', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep2), 'Hep2', merged_samples$Hep_clusters) 
merged_samples$Hep_clusters = ifelse(merged_samples$Hep_clusters %in% as.character(Hep3), 'Hep3', merged_samples$Hep_clusters) 
table(merged_samples$Hep_clusters)
####################################################

#Heat map 1  - human only
CV_list = c("Akr1c1", "Cyp7a1", "Cyp27a1", "Ahr", "Fasn", "Glul", "Notum", "Rcan1", 'Por')
#PP_list = c("Alb", "Apoa1", "Apoc3", "Apoe", "Cyp2f4", "Ass1", "Acly", "Hal", "Sds", "Orm1")
PP_list = c("Alb", "Apoa1", "Apoc1", "Apoc3", "Apoe", "Cps1", "Fabp1", 
            "Orm1", "Pigr", "Serpina1", "Tf", "Ttr")
markers_list = c(CV_list, PP_list)
####################################################
#Heat map 2  - mix
CV_list =  c("Akr1c1", "Cyp7a1", "Cyp27a1", "Ahr", "Slco1b2", "Fasn", "Glul", "Notum", "Rcan1", "Cyp2e1")
PP_list =  c("Alb", "Apoa1", "Apoc1", "Apoc3", "Apoe", "Tf", "Ttr", 
             "Orm1", "Arg1", "Pck1", "Sds", "Cyp2f4", "Hal")
markers_list = c(CV_list, PP_list)
####################################################
# heatmap 3 - mouse only
CV_list =  c("Axin2", "Cyp27a1", "Cyp2e1", "Cyp7a1", "Fasn", "Glul", "Slco1b2", "Psmd4")
PP_list =   c("Arg1", "Pck1", "Sds", "Cyp2f4", "Hal", "Acly", "Ass1", "Shdh", "Srebf1", "Acly")
markers_list = c(CV_list, PP_list)

################################

CV_list=c("Akr1c1", "Cyp7a1", "Cyp27a1", "Ahr", "Glul", "Notum", "Rcan1", "Por", "Slco1b2", "Cyp1a2", "Cyp2e1")
IZ_list=c("Saa4", "Hint1", "Cyp8b1", "Cyp2e1", "Cox7c", "Fabp1")
PP_list=c("Alb", "Apoa1", "Apoc1", "Apoc3", "Apoe", "Serpina1", "Ttr", "Orm1", "Acly", "Ass1", "Cyp2f4", "Sds", "Hal")
markers_list = c(CV_list, IZ_list, PP_list)

################################

PP1_list= c("Sds", "Hal", "Cyp2f4", "Arg1", "Acly", "Ass1", "Gls2", "Agxt", "Uroc1", "Gldc", "Gls2")
PP2_list= c("Apoc3", "Serpina1", "Apoc1", "Apoe", "Itih4", "Apoa1", "Ttr", "Tf", "Alb", "Pigr", "Orm1", "Rpl3", "Fads1", "Aldh1b1", "Srd5a1", "Hsd17b13")
IZ_list= c("Saa4", "Hint1", "Cyp8b1", "Cyp2e1", "Cox7c", "Fabp1")
CV_list= c("Notum", "Cyp27a1", "Fabp7", "Akr1c1", "Gsta5", "Slc22a1", "Aox3",
           "Sult1e1", "Fmo1", "Oat", "Ahr", "Cyp7a1", "Glul", "Rhbg", "Cyp2e1", "Cyp1a2", "Por")
markers_list = c(CV_list, IZ_list, PP2_list, PP1_list)

##########################
#### markers based on Spatial data PCA factors
PC_biased = c("Scd", "Ca3", "Ces1d", "Ugt1a5", "Gulo", "Fetub", "Gstt3", "Cyp2c22", 
              "Xpnpep2", "Cdh17", "Cyp1a2", "Cyp2c11", "C1r", "Cyp2e1", "AABR07047899.1",
              "AABR07048463.1", "AABR07048474.1", "LOC259244", "LOC100910877", "LOC100360095", 
              "Slc27a5", "Gstm1", "Cyp7a1", "Cela1", "Slc13a3", "Dhrs7l1", 
              "Mup4", "Akr7a3", "Oat", "Sord", "Rgn", "Ahr", "Fmo1", "Avpr1a", "Slc22a1", 
              "Sult1e1", "Csad", "Nr1i3", "Akr1c2", "Cyp27a1", "Glul")


PP_biased = c("Apoa4", "Orm1", "Hmgcs1", "Cdh1", "Apoc2", "Cyp4a2", "Atp1a1", "Wfdc2", "Srd5a1",
             "Spink1l", "Itih4", "Wfdc21", "Gjb2", "Apoc3", "Gpx1", "Apoa1", "Alb", "Rarres2", 
             "Cfi", "Fabp1", "Ctsh", "Uroc1", "Etnppl", "Hsd17b13", "Slc38a4", "Sult2a6", "Ass1",
             "Gldc", "Arg1", "Agxt", "Cyp2c7", "Crot", "Ftcd", "Sds", "Hal", "Mfsd2a", "Cyp4a2", 
             "Cps1", "Cyp2c7", "Tdo2", "Slc25a47", "Gls2", "Bhmt")

markers_list = c(PC_biased,PP_biased)

################################
pericentral = c("Glul","Fmo1", "Oat", "Aox3", "Akr1c1", "Akr1c2", "Csad", "Sult1e1", "Fabp7", "Mup4", "Slc27a5", 
                "Slc22a1", "Gsta5", "Cyp2e1", "Slc13a3", "Gulo", "Cyp1a2", "Scd", "Ces1d", "Ugt1a5", "Cyp2c11", "Cyp2c22")
non_monotonic = c("Por", "Notum", "Cyp27a1", "Cox7c", "Cyp8b1", "Hamp")
periportal = c("Pigr", "Itih4", "Apoc3", 
                "Serpina1", "Fads1", "Rpl3", "Cfi", "Gpx1", "Fabp1", "Apoa1", "Srd5a1", "Spink1l", "Orm1", "Tf", 
                "Alb", "Apoa4", "Gjb2", "Aldh1b1", "Cyp4a2", "Acly", "Hmgcs1", "Sult2a6", "Cyp2c7", "Slc38a4",
                "Ass1", "Cps1", "Hsd17b13", "Etnppl", "Ftcd", "Sds", "Mfsd2a", "Sult2a6", "Gldc", "Arg1", 
                "Hal", "Uroc1", "Agxt", "Gldc", "Gls2")
markers_list = c(pericentral, non_monotonic, periportal)
################################
##### Generating average expression based heatmap

get_scaled_by_gene <- function(x){
  y <- scale(t(x), center = TRUE, scale = TRUE)
  y[is.nan(y)] <- 0
  y2 <- t(y)
  y2 <- as.data.frame(y2)
  return(y2)
}


merged_samples_subset = merged_samples[, merged_samples$Hep_clusters %in% c('Hep0', 'Hep1', 'Hep2', 'Hep3')]

#merged_samples_subset2 <- CreateAssayObject(GetAssayData(merged_samples_subset)[rownames(merged_samples_subset) %in% markers_list,])
merged_samples_subset <- merged_samples_subset[rownames(merged_samples_subset) %in% markers_list,]

dim(merged_samples_subset)
#merged_samples_subset2$manual_annotation = merged_samples_subset$manual_annotation
cluster_names = as.character(merged_samples_subset$SCT_snn_res.2.5 )
cluster_names_types = names(table(cluster_names))

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  GetAssayData(merged_samples_subset, 'data')[,cluster_names == a_cluster_name] 
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)


## calculate the average expression of each gene in each cluster
cluster_average_exp <- lapply(cluster_expression, function(x){
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
lapply(cluster_average_exp, dim)

## Concatenate all the clusters together to make a matrix
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
colnames(cluster_average_exp_df) = paste0(names(cluster_average_exp))
head(cluster_average_exp_df)


## scale and center all the genes in the matrix
cluster_average_exp_df_scaled <- get_scaled_by_gene(cluster_average_exp_df) ## scaling over the clusters of interest
colnames(cluster_average_exp_df_scaled)
rownames(cluster_average_exp_df_scaled)

anno_df = c(rep(NA, ncol(cluster_average_exp_df_scaled)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Hep0) anno_df[i] <<- 'Hep0 (PP-like 1)'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Hep1) anno_df[i] <<- 'Hep1 (PP-like 2)'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Hep2) anno_df[i] <<- 'Hep2 (CV-like 1)'
  if (colnames(cluster_average_exp_df_scaled)[i] %in% Hep3) anno_df[i] <<- 'Hep3 (CV-like 2)'

})


anno_df_genes = c(rep(NA, nrow(cluster_average_exp_df_scaled)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (rownames(cluster_average_exp_df_scaled)[i] %in% CV_list) anno_df_genes[i] <<- 'CV'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% IZ_list) anno_df_genes[i] <<- 'IZ'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% PP2_list) anno_df_genes[i] <<- 'PP2'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% PP1_list) anno_df_genes[i] <<- 'PP1'
})

anno_df_genes = c(rep(NA, nrow(cluster_average_exp_df_scaled)))
sapply(1:ncol(merged_samples_subset) , function(i){
  if (rownames(cluster_average_exp_df_scaled)[i] %in% pericentral) anno_df_genes[i] <<- 'Pericentrally-enriched'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% non_monotonic) anno_df_genes[i] <<- 'Non-monotonic'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% periportal) anno_df_genes[i] <<- 'Periportally-enriched'
})

sapply(1:ncol(merged_samples_subset) , function(i){
  if (rownames(cluster_average_exp_df_scaled)[i] %in% pericentral) anno_df_genes[i] <<- 'Pericentrally-enriched'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% non_monotonic) anno_df_genes[i] <<- 'Non-monotonic'
  if (rownames(cluster_average_exp_df_scaled)[i] %in% periportal) anno_df_genes[i] <<- 'Periportally-enriched'
})
length(anno_df_genes)

##### ordering clusters
cluster_orders = as.character(c(Hep3, Hep2, Hep0, Hep1))
anno_df = data.frame(cluster=colnames(cluster_average_exp_df_scaled),annot=anno_df)

cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[,cluster_orders]
anno_df = anno_df[match(colnames(cluster_average_exp_df_scaled), anno_df$cluster),]
head(anno_df)
dim(cluster_average_exp_df_scaled)



#### ordering genes
gene_orders = c(CV_list, IZ_list, PP2_list, PP1_list) 
gene_orders = c(PC_biased, PP_biased)
gene_orders = markers_list

gene_orders = unique(gene_orders[gene_orders %in% rownames(cluster_average_exp_df_scaled)])
anno_df_genes = data.frame(genes=row.names(cluster_average_exp_df_scaled),annot=anno_df_genes)


cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[gene_orders,]
cluster_average_exp_df_scaled = cluster_average_exp_df_scaled[,cluster_orders]


dim(cluster_average_exp_df_scaled)
anno_df_genes = anno_df_genes[match(row.names(cluster_average_exp_df_scaled), anno_df_genes$genes),]



dim(anno_df)
dim(anno_df_genes)

row.names(anno_df) = anno_df$cluster
row.names(anno_df_genes) = anno_df_genes$genes

head(anno_df)
head(anno_df_genes)

anno_df = anno_df[,-1,drop=F]
anno_df_genes = anno_df_genes[,-1,drop=F ]

dim(cluster_average_exp_df_scaled)


rownames(cluster_average_exp_df_scaled) %in% rownames(anno_df_genes)
colnames(cluster_average_exp_df_scaled) %in% rownames(anno_df)



cluster_average_exp_df_scaled = t(cluster_average_exp_df_scaled)

colnames(anno_df) = 'Hep.group'
colnames(anno_df_genes) = 'Zonation'


hep_colors=c('royalblue3','steelblue2' ,'tomato1', 'red2')
names(hep_colors) = unique(anno_df$Hep.group)

zonation_colors=c('skyblue3','orange' ,'red3')
names(zonation_colors) = unique(anno_df_genes$Zonation)

pheatmap(cluster_average_exp_df_scaled, cluster_rows = FALSE,
         cluster_cols = FALSE, 
         annotation_row = data.frame(anno_df),
         annotation_col = data.frame(anno_df_genes),
         color= inferno(20, direction = +1),  
         #border_color= TRUE, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         annotation_colors = list(Hep.group=hep_colors, Zonation=zonation_colors),
         annotation_legend = TRUE
)




