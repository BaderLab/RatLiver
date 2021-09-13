source('Codes/Functions.R')
Initialize()

################################
##### merging all the information related to the set-1
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

top_n = 30
get_head_genes <- function(df){df$genes[1:top_n] }
get_tail_genes <- function(df){ df$genes[(nrow(df)-top_n):nrow(df)] }

get_head_10genes <- function(df){df$genes[1:10] }
get_tail_10genes <- function(df){ df$genes[(nrow(df)-10):nrow(df)] }


old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
set1_data <- your_scRNAseq_data_object
set1_data@meta.data$cluster = as.character(sCVdata_list$res.0.6@Clusters)

varimax_res_set1 <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
score_set1 <- data.frame(varimax_res_set1$rotScores)
colnames(score_set1) = paste0('Varimax_', 1:ncol(score_set1))
loading_set1 <- varimax_res_set1$rotLoadings

set1_info <- read.csv('figure_panel/set-1-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_data@meta.data$umi = colnames(set1_data)
meta.data.label = merge(set1_data@meta.data, set1_info, by.x='cluster', by.y='clusters')
meta.data.label <- meta.data.label[match(colnames(set1_data),meta.data.label$umi),]
head(meta.data.label)

set1_df <- data.frame(umi=rownames(score_set1), 
                 strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2)), 
                 cluster = set1_data$cluster, 
                 meta.data.label,
                 score_set1)
head(set1_df)

varimax_5_set1 = data.frame(genes=rownames(loading_set1),load=loading_set1[,5])
varimax_5_set1 = varimax_5_set1[order(varimax_5_set1$load, decreasing = T),]
head(varimax_5_set1)

varimax_15_set1 = data.frame(genes=rownames(loading_set1),load=loading_set1[,15])
varimax_15_set1 = varimax_15_set1[order(varimax_15_set1$load, decreasing = T),]
head(varimax_15_set1)

top_n = 30
get_head_genes <- function(df){df$genes[1:top_n] }
get_tail_genes <- function(df){ df$genes[(nrow(df)-top_n):nrow(df)] }


## on the DA side: RT1-A2
get_head_genes(varimax_5_set1)[get_head_genes(varimax_5_set1) %in% get_tail_genes(varimax_15_set1)]
### 2 lewis genes: "Tmem176a" "RT1-A1" 
get_tail_genes(varimax_5_set1)[get_tail_genes(varimax_5_set1) %in% get_head_genes(varimax_15_set1)]
 

################################
##### merging all the information related to the set-2

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_labelCor.RData"
load(new_data_scCLustViz_object)
set2_data <- your_scRNAseq_data_object
set2_data@meta.data$cluster = as.character(sCVdata_list$res.0.6@Clusters)

varimax_res_set2 <- readRDS('Results/new_samples/varimax_rotated_object_new.rds')
score_set2 <- data.frame(varimax_res_set2$rotScores)
colnames(score_set2) = paste0('Varimax_', 1:ncol(score_set2))
loading_set2 <- varimax_res_set2$rotLoadings

set2_info <- read.csv('figure_panel/set-2-final-info.csv')
colnames(set2_info)[1] = 'clusters'
set2_data@meta.data$umi = colnames(set2_data)
meta.data.label = merge(set2_data@meta.data, set2_info, by.x='cluster', by.y='clusters')
meta.data.label <- meta.data.label[match(colnames(set2_data),meta.data.label$umi),]
head(meta.data.label)
set2_data@meta.data = meta.data.label

set2_df <- data.frame(umi=rownames(score_set2), 
                      strain = unlist(lapply(str_split(colnames(set2_data),'_'), '[[', 2)), 
                      cluster = set2_data$cluster, 
                      meta.data.label,
                      score_set2)
head(set2_df)
set2_data@meta.data = set2_df



varimax_5_set2 = data.frame(genes=rownames(loading_set2),load=loading_set2[,5])
varimax_5_set2 = varimax_5_set2[order(varimax_5_set2$load, decreasing = T),]
head(varimax_5_set2)
tail(varimax_5_set2)

varimax_16_set2 = data.frame(genes=rownames(loading_set2),load=loading_set2[,16])
varimax_16_set2 = varimax_16_set2[order(varimax_16_set2$load, decreasing = T),]
head(varimax_16_set2)



## Apoc1 
get_head_genes(varimax_5_set2)[get_head_genes(varimax_5_set2) %in% get_head_genes(varimax_5_set1)]
## "RT1-CE10" "Tsc22d3"  "RT1-A2"  
get_tail_genes(varimax_5_set2)[get_tail_genes(varimax_5_set2) %in% get_head_genes(varimax_5_set1)]

## Cth
get_head_genes(varimax_16_set2)[get_head_genes(varimax_16_set2) %in% get_head_genes(varimax_5_set1)]
## Tpt1
get_head_genes(varimax_16_set2)[get_head_genes(varimax_16_set2) %in% get_tail_genes(varimax_5_set1)]
## "Apoa4"  "Cyp2d5"
get_tail_genes(varimax_16_set2)[get_tail_genes(varimax_16_set2) %in% get_tail_genes(varimax_5_set1)]

################################




################################
##### merging all the information related to the set-2

new_immune_scCLustViz_object <- 'Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub_c17Included_labelCor.RData'
load(new_immune_scCLustViz_object)
set2_immune_data = your_scRNAseq_data_object

varimax_res_immune = readRDS('Results/new_samples/immune_varimax_results_cl17Inc.rds')
score_set2_immune <- data.frame(varimax_res_immune$rotScores)
rownames(score_set2_immune) = colnames(set2_immune_data)
colnames(score_set2_immune) = paste0('Varimax_', 1:ncol(score_set2_immune))
loading_set2_immune <- varimax_res_immune$rotLoadings

set2_immune_info <- read.csv('figure_panel/set-2-immunesub-final-info.csv')
colnames(set2_immune_info)[1] = 'clusters'
set2_immune_data@meta.data$umi = colnames(set2_immune_data)
meta.data.label = merge(set2_immune_data@meta.data, set2_immune_info, by.x='cluster', by.y='clusters')
meta.data.label <- meta.data.label[match(colnames(set2_immune_data),meta.data.label$umi),]
head(meta.data.label)

set2_immune_df <- data.frame(umi=colnames(set2_immune_data), 
                      strain = set2_immune_data$strain,
                      cluster = set2_immune_data$cluster, 
                      meta.data.label,
                      score_set2_immune)
head(set2_immune_df)



varimax_2_immune = data.frame(genes=rownames(loading_set2_immune),load=loading_set2_immune[,2])
varimax_2_immune = varimax_2_immune[order(varimax_2_immune$load, decreasing = T),]
head(varimax_2_immune)

ggplot(set2_immune_df,aes(Varimax_1,Varimax_2,color=strain))+geom_point()+theme_classic()


## Apoc1
get_head_genes(varimax_2_immune)[get_head_genes(varimax_2_immune) %in% get_head_genes(varimax_5_set1)]
## "Ttr , Apoc3, AABR07048474.1, LOC680406, Rup2", "Fabp1", "Apoc1", "Serpina1", "Rbp4"  
## Alb, Hp, Ugt2b, Gc, LOC297568, Apoa2, Fgb, Apoe, Tf, LOC259244, Aldob   
get_head_genes(varimax_2_immune)[get_head_genes(varimax_2_immune) %in% get_head_genes(varimax_5_set2)]
## "Fabp1"        "Serpina1"     "Ugt2b"        "LOC100360791" "Cyp2c6v1"    
get_head_genes(varimax_2_immune)[get_head_genes(varimax_2_immune) %in% get_tail_genes(varimax_5_set1)]
## "Ambp"  "Hpx"   "Fetub"
get_head_genes(varimax_2_immune)[get_head_genes(varimax_2_immune) %in% get_tail_genes(varimax_16_set2)]

# "RT1-A2"
get_tail_genes(varimax_2_immune)[get_tail_genes(varimax_2_immune) %in% get_head_genes(varimax_5_set1)]
## "RGD1559482" "Rgs1"       "RGD1562667" "RT1-A2"    
get_tail_genes(varimax_2_immune)[get_tail_genes(varimax_2_immune) %in% get_tail_genes(varimax_15_set1)]
## "Car2"       "Oasl2"      "RGD1562667" "RT1-A2"
get_tail_genes(varimax_2_immune)[get_tail_genes(varimax_2_immune) %in% get_tail_genes(varimax_5_set2)]

################################



strain_genes = c('RT1-A2', 'Tmem176a', 'RT1-A1', 'Apoc1', "RT1-CE10", "Tsc22d3", "RT1-A2", 'Cth', 'Tpt1', "Apoa4" , "Cyp2d5",
  "Ttr" , 'Apoc3', 'AABR07048474.1', 'LOC680406', "Rup2", "Fabp1", "Apoc1", "Serpina1", "Rbp4",
  'Alb', 'Hp', 'Ugt2b', 'Gc', 'LOC297568', 'Apoa2', 'Fgb', 'Apoe', 'Tf', 'LOC259244', 'Aldob',
  "Fabp1", "Serpina1", "Ugt2b", "LOC100360791", "Cyp2c6v1",
  "Ambp", "Hpx", "Fetub", "RGD1559482", "Rgs1", "RGD1562667", "RT1-A2", 
  "Car2", "Oasl2", "RGD1562667", "RT1-A2")

strain_genes = unique(strain_genes)
write.csv(data.frame(strain_genes),'Results/strain_specific_genes.csv')


strain_genes = c(get_tail_genes(varimax_15_set1), 
        get_tail_genes(varimax_5_set1), 
        get_head_genes(varimax_15_set1), 
        get_head_genes(varimax_5_set1))

length(unique(strain_genes))
length(strain_genes)
#strain_genes = c(strain_genes, 'Ptprc', 'Actb')
### visualize them as a cell-type specific boxplot

###### a quick of the validated genes 
#strain_genes = c("Itgal", "AC114233.2", "Timp2", "RT1-A1" ) ### LEW genes
strain_genes = x
#strain_genes = c("RGD1562667", "Siglec5", "RGD1559482", "Ly6al" ) ## DA genes
strain_genes = y
######

geneExpData_set1 = GetAssayData(set1_data[strain_genes,])
geneExpData_set1 = t(geneExpData_set1)

geneExpData_set2 = GetAssayData(set2_data[strain_genes,])
geneExpData_set2 = t(geneExpData_set2)


#pdf('Plots/Expression_boxplot_strainSpecificGenes.pdf',width = 13,height = 7.5)
pdf('Plots/Expression_boxplot_validatedWithDE_DA_2.pdf',width = 13,height = 7.5)
for(i in 1:length(strain_genes)){
  gene_name = colnames(geneExpData_set1)[i]
  
  df_set1 = data.frame(gene_set1=geneExpData_set1[,i],
                       label=set1_df$label,
                       strain=set1_df$strain)
  
  df_set1$label = gsub('Inflammatory','Inf',df_set1$label)
  df_set1$label = gsub('and','&',df_set1$label)

  p1=ggplot(df_set1, aes(x=label, y=gene_set1, fill=strain))+geom_boxplot()+theme_classic()+
    scale_fill_manual(values = c("#CC79A7","#0072B2"))+
    theme_classic()+ylab(gene_name)+ggtitle('set1')+
    theme(text = element_text(size=13),
          axis.text.x = element_text(size=10,angle=90,color='black'),
          legend.title = element_blank())+xlab('') # legend.position = "none"
  
  p2=ggplot(df_set1, aes(x=strain, y=gene_set1, fill=strain))+geom_boxplot()+theme_classic()+
    scale_fill_manual(values = c("#CC79A7","#0072B2"))+
    theme_classic()+ylab(gene_name)+ggtitle('set1')+
    stat_summary(fun.data=data_summary)+stat_compare_means()+#method = "t.test"
    theme(text = element_text(size=13),
          axis.text.x = element_text(size=10,angle=90,color='black'),
          legend.title = element_blank())+xlab('') # legend.position = "none"
  
  
  
  ##################################
  gene_name = colnames(geneExpData_set2)[i]
  
  df_set2 = data.frame(gene_set2=geneExpData_set2[,i],
                       label=set2_df$label,
                       strain=set2_df$strain)
  
  df_set2$label = gsub('Inflammatory','Inf',df_set2$label)
  table(df_set2$label)
  
  p3=ggplot(df_set2, aes(x=label, y=gene_set2, fill=strain))+geom_boxplot()+theme_classic()+
    scale_fill_manual(values = c("#CC79A7","#0072B2"))+
    theme_classic()+ylab(gene_name)+ggtitle('set2')+
    theme(text = element_text(size=13),
          axis.text.x = element_text(size=10,angle=90,color='black'),
          legend.title = element_blank())+xlab('') # legend.position = "none"
  
  p4=ggplot(df_set2, aes(x=strain, y=gene_set2, fill=strain))+geom_boxplot()+theme_classic()+
    scale_fill_manual(values = c("#CC79A7","#0072B2"))+
    theme_classic()+ylab(gene_name)+ggtitle('set2')+
    stat_summary(fun.data=data_summary)+stat_compare_means()+#method = "t.test"
    theme(text = element_text(size=13),
          axis.text.x = element_text(size=10,angle=90,color='black'),
          legend.title = element_blank())+xlab('') # legend.position = "none"
  
  ##################################
  gridExtra::grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
  
  
}
dev.off()



####################################################################
############ Running DE analysis between the strains   
####################################################################
cell_types_set2 = names(table(set2_data$label))
set2_data_split = sapply(1:length(cell_types_set2), function(i) set2_data[,set2_data$label==cell_types_set2[i]], simplify = F)
names(set2_data_split) = cell_types_set2
lapply(set2_data_split, dim)
sapply(1:length(set2_data_split), function(i) Idents(set2_data_split[[i]]) <<- set2_data_split[[i]]$strain)

strain_vector = names(table(set2_data$strain))
DAvsLEW_cellTypeMarkers = sapply(1:length(set2_data_split), function(i)  FindMarkers(set2_data_split[[i]], 
                                                ident.1='DA', #DA
                                                ident.2='LEW', #LEW
                                                logfc.threshold = 0,
                                                min.pct=0,
                                                min.cells.feature = 1,
                                                min.cells.group = 1), simplify = FALSE)

names(DAvsLEW_cellTypeMarkers) = names(set2_data_split)

DAvsLEW_cellTypeMarkers = readRDS('Results/DAvsLEW_cellTypeMarkers_set2.rds')
lapply(DAvsLEW_cellTypeMarkers, head)

DAvsLEW_cellTypeMarkers.ord = lapply(DAvsLEW_cellTypeMarkers, function(x) {x=x[order(x$avg_log2FC,decreasing = T),]; x=x[x$p_val_adj<0.01,];x})
toInclude = unlist(lapply(DAvsLEW_cellTypeMarkers.ord, function(x) nrow(x)))>5
DAvsLEW_cellTypeMarkers.ord =DAvsLEW_cellTypeMarkers.ord[toInclude]

top_DE_num = 20
DAvsLEW_genes = lapply(DAvsLEW_cellTypeMarkers.ord, function(x) rownames(x)[c(1:top_DE_num)])
DAvsLEW_genes.df = data.frame(table(unlist(DAvsLEW_genes)))
DAvsLEW_genes.df[order(DAvsLEW_genes.df$Freq,decreasing = T),]

LEWvsDA_genes = lapply(DAvsLEW_cellTypeMarkers.ord, function(x) rownames(x)[(nrow(x)-top_DE_num):nrow(x)])
LEWvsDA_genes.df = data.frame(table(unlist(LEWvsDA_genes)))
LEWvsDA_genes.df[order(LEWvsDA_genes.df$Freq,decreasing = T),]


Reduce(intersect, list(strain_genes,DAvsLEW_genes.df$Var1)) ### DA specific genes
Reduce(intersect, list(strain_genes,LEWvsDA_genes.df$Var1)) ### LEW specific genes 

markersToCheck = c('Macro', 'Cd163', 'Hmox1', 'Vsig4', 'Cd5l', 'Mrc1')
markersToCheck %in% DAvsLEW_genes.df$Var1
markersToCheck[markersToCheck %in% LEWvsDA_genes.df$Var1]

'Mrc1' %in% LEWvsDA_genes.df$Var1
'Mrc1' %in% DAvsLEW_genes.df$Var1
####################################################################
### maybe run a generic (non-cell type specific version as well)
Idents(set2_data) = set2_data$strain
generalMarkers_DAvsLEW = FindMarkers(set2_data, 
            ident.1='DA', #DA
            ident.2='LEW', #LEW
            logfc.threshold = 0,
            min.pct=0,
            min.cells.feature = 1,
            min.cells.group = 1)

generalMarkers_DAvsLEW = readRDS('Results/generalMarkers_DAvsLEW_set2.rds')

generalMarkers_DAvsLEW = generalMarkers_DAvsLEW[order(generalMarkers_DAvsLEW$avg_log2FC,decreasing = T),]
generalMarkers_DAvsLEW = generalMarkers_DAvsLEW[generalMarkers_DAvsLEW$p_val_adj<0.01,]

rownames(generalMarkers_DAvsLEW)[1:10]
numRow = nrow(generalMarkers_DAvsLEW)
rownames(generalMarkers_DAvsLEW)[(numRow-10):numRow]


#######################################
############# checking the enrichment of set1 varimax 5 and 15 in the set2 cells using Ucell enrrrichment
markers <- list()
markers$var15head = get_head_10genes(varimax_15_set1)
markers$var15tail = get_tail_10genes(varimax_15_set1)

markers$var5head = get_head_10genes(varimax_5_set1)
markers$var5tail = get_tail_10genes(varimax_5_set1)

lapply(markers, head)
my.matrix <- GetAssayData(set2_data)

### calculating enrichment score for each signature
scores <- data.frame(ScoreSignatures_UCell(my.matrix, features=markers))
rownames(scores) == colnames(set2_data)

set2.umap.df = data.frame(Embeddings(set2_data, 'umap'), 
                          scores,
                          strain=set2_data$strain,
                          label=set2_data$label)
head(set2.umap.df)

##### projecting varimax 15 on the set2 map
ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var15head_UCell))+geom_point(alpha=0.7,size=1)+
  theme_classic()+scale_color_viridis(option = 'magma',direction = -1)
ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var15tail_UCell))+geom_point(alpha=0.7,size=1)+
  theme_classic()+scale_color_viridis(option = 'magma',direction = -1)

##### projecting varimax 5 on the set2 map
ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var5head_UCell))+geom_point(alpha=0.7,size=1)+
  theme_classic()+scale_color_viridis(option = 'magma',direction = -1)
ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var5tail_UCell))+geom_point(alpha=0.7,size=1)+
  theme_classic()+scale_color_viridis(option = 'magma',direction = -1)


### is varimax 15 strain specific
ggplot(set2.umap.df, aes( x=strain,y=var15head_UCell,fill=strain))+geom_boxplot()+theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))
ggplot(set2.umap.df, aes( x=strain,y=var15tail_UCell,fill=strain))+geom_boxplot()+theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))
### is varimax 5 strain specific
ggplot(set2.umap.df, aes( x=strain,y=var5head_UCell,fill=strain))+geom_boxplot()+theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))
ggplot(set2.umap.df, aes( x=strain,y=var5tail_UCell,fill=strain))+geom_boxplot()+theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))

table(set2.umap.df$label[set2.umap.df$var5head_UCell>0.3])
#### checking if varimax-5 is valid in the Hep population of the set2 map
set2.umap.df.sub = set2.umap.df[str_sub(set2.umap.df$label, 1,3)=='Hep',]
ggplot(set2.umap.df.sub, aes( x=strain,y=var5head_UCell,fill=strain))+geom_boxplot()+theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))
ggplot(set2.umap.df.sub, aes( x=strain,y=var5tail_UCell,fill=strain))+geom_boxplot()+theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))


###### checking consistency with the DE genes:
get_head_10genes(varimax_15_set1)[get_head_10genes(varimax_15_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes
get_tail_10genes(varimax_15_set1)[get_tail_10genes(varimax_15_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes

get_head_10genes(varimax_5_set1)[get_head_10genes(varimax_5_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes
get_tail_10genes(varimax_5_set1)[get_tail_10genes(varimax_5_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes

get_head_10genes(varimax_5_set1)[get_head_10genes(varimax_5_set1) %in% LEWvsDA_genes.df$Var1] ### should not give any results
get_tail_10genes(varimax_5_set1)[get_tail_10genes(varimax_5_set1) %in% DAvsLEW_genes.df$Var1] ### should not give any results

x=get_head_genes(varimax_15_set1)[get_head_genes(varimax_15_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes
y=get_tail_genes(varimax_15_set1)[get_tail_genes(varimax_15_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes
final_strain_markers = c(x,y)
final_strain_markers = final_strain_markers[!final_strain_markers %in% 'RT1-CE10']


x=get_head_genes(varimax_5_set1)[get_head_genes(varimax_5_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes
y=get_tail_genes(varimax_5_set1)[get_tail_genes(varimax_5_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes
final_strain_markers = c(y,x)
final_strain_markers = final_strain_markers[!final_strain_markers %in% c('Apoc1', 'Scd')]


final_strain_markers =  c("Itgal","AC114233.2","Timp2","RT1-A1","Tmem176a" ,
                          "RT1-A2" ,"RGD1562667","Siglec5","RGD1559482","Ly6al")

DoHeatmap((set2_data), features = final_strain_markers,slot='counts')
DotPlot(set2_data,features =final_strain_markers) + 
  RotatedAxis()+ xlab('Markers')+ylab('')+theme(axis.text.x = element_text(size=11))

table(set2_data$label)were
set2_data_tmp = set2_data[,set2_data$cluster %in% as.character(c(17,10,11,8,12,7,13))] # all immune
set2_data_tmp = set2_data[,set2_data$cluster %in% as.character(c(11, 8, 13))] ## Mac
set2_data_tmp = set2_data[,set2_data$cluster %in% as.character(c(0, 1, 2, 3, 9, 14, 15))] ## Hep

DotPlot(set2_data_tmp,features =final_strain_markers) + 
  RotatedAxis()+ xlab('Markers')+ylab('')+theme(axis.text.x = element_text(size=11))


Idents(set2_data_tmp) = set2_data_tmp$label

DotPlot(set2_data_tmp, split.by = "strain",
        features = final_strain_markers) + 
  RotatedAxis()+ xlab('Markers')+ylab('')+
  theme(axis.text.x = element_text(size=11))
#split.by = "groups"

