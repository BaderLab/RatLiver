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
set1_data@meta.data$strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2))

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
set1_data@meta.data = meta.data.label
table(set1_data$label)
set1_data$sample = unlist(lapply(str_split(set1_data$umi, '_'), function(x) paste0(x[2:3], collapse = '_')))
  
set1_df <- data.frame(umi=rownames(score_set1), 
                 strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2)), 
                 cluster = set1_data$cluster, 
                 meta.data.label,
                 score_set1)
head(set1_df)

###### Varimax 5 
varimax_5_set1 = data.frame(genes=rownames(loading_set1),load=loading_set1[,5])
varimax_5_set1 = varimax_5_set1[order(varimax_5_set1$load, decreasing = T),]
rownames(varimax_5_set1) = NULL

varimax_5_set1.vis = varimax_5_set1
colnames(varimax_5_set1.vis) = c('Genes', 'Loading score')
varimax_5_set1.vis$`Loading score` = round(varimax_5_set1$load, 3)
gridExtra::grid.table(head(varimax_5_set1.vis,10))
dev.off()
varimax_5_set1.vis = varimax_5_set1.vis[order(varimax_5_set1$load, decreasing = F),]
gridExtra::grid.table(head(varimax_5_set1.vis,10))

###### Varimax 15 
varimax_15_set1 = data.frame(genes=rownames(loading_set1),load=loading_set1[,15])
varimax_15_set1 = varimax_15_set1[order(varimax_15_set1$load, decreasing = T),]
rownames(varimax_15_set1) = NULL

varimax_15_set1.vis = varimax_15_set1
colnames(varimax_15_set1.vis) = c('Genes', 'Loading score')
varimax_15_set1.vis$`Loading score` = round(varimax_15_set1$load, 3)
gridExtra::grid.table(head(varimax_15_set1.vis,10))
dev.off()
varimax_15_set1.vis.rev = varimax_15_set1.vis[order(varimax_15_set1$load, decreasing = F),]
gridExtra::grid.table(head(varimax_15_set1.vis.rev,10))

varimax_15_set1.Lew = varimax_15_set1.vis[1:25,]
varimax_15_set1.Lew$strain = 'LEW'
varimax_15_set1.DA = varimax_15_set1.vis.rev[1:2000,]
varimax_15_set1.DA$strain = 'DA'
write.csv(rbind(varimax_15_set1.Lew, varimax_15_set1.DA), 'figure_panel/set1-varimax15-topGenes.csv')
# scp delaram@192.168.142.161:~/RatLiver/figure_panel/set1-varimax15-topGenes.csv ~/Desktop/

## on the DA side: RT1-A2
get_head_genes(varimax_5_set1)[get_head_genes(varimax_5_set1) %in% get_tail_genes(varimax_15_set1)]
### 2 lewis genes: "Tmem176a" "RT1-A1" 
get_tail_genes(varimax_5_set1)[get_tail_genes(varimax_5_set1) %in% get_head_genes(varimax_15_set1)]
 
dev.off()


#### checking the expression pattern of top loaded markers over UMAP
i = nrow(varimax_15_set1) -4
gene2check = varimax_15_set1$genes[i]
gene2check
umap.df = data.frame(getEmb(set1_data, 'umap'),gene=GetAssayData(set1_data)[gene2check,])
ggplot(umap.df, aes(UMAP_1,UMAP_2, color=gene))+geom_point(alpha=0.1)+theme_classic()+scale_color_viridis(direction = -1)
set1_df$gene = umap.df$gene
ggplot(set1_df, aes(Varimax_1, Varimax_15, color=gene))+geom_point(size=1.2,alpha=0.6)+
  theme_classic()+scale_color_viridis(direction = -1)+ggtitle(gene2check)+
  theme(text = element_text(size=17.5),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())
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




########## various ways to select the genes of interest
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





#####################################################################################################################
############# checking the enrichment of set1 varimax 5 and 15 in the set2 cells using Ucell enrrrichment #######
#####################################################################################################################
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
ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var15head_UCell))+geom_point(alpha=0.7,size=1.3)+
  theme_classic()+scale_color_viridis('Lew macrophage-\nspecific signature score\n(based on set-1 Varimax-15)',option = 'magma',direction = -1)+
  ggtitle('Enrichment score of Var-15 LEW macrophage-specific signatures \n(based on set1 samples) over set2 samples UMAP')

ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var15tail_UCell))+geom_point(alpha=0.7,size=1.3)+
  theme_classic()+scale_color_viridis('DA macrophage-\nspecific signature score\n(based on set-1 Varimax-15)',option = 'magma',direction = -1)+
  ggtitle('Enrichment score of Var-15 DA macrophage-specific signatures \n(based on set1 samples) over set2 samples UMAP')

##### projecting varimax 5 on the set2 map
ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var5head_UCell))+geom_point(alpha=0.7,size=1.3)+
  theme_classic()+scale_color_viridis('DA hepatocyte-\nspecific signature score\n(based on set-1 Varimax-5)',option = 'magma',direction = -1)+
  ggtitle('Enrichment score of Var-5 DA hepatocyte-specific signatures \n(based on set1 samples) over set2 samples UMAP')

ggplot(set2.umap.df, aes( x=UMAP_1,y=UMAP_2,color=var5tail_UCell))+geom_point(alpha=0.7,size=1.3)+
  theme_classic()+scale_color_viridis('Lew hepatocyte-\nspecific signature score\n(based on set-1 Varimax-5)',option = 'magma',direction = -1)+
  ggtitle('Enrichment score of Var-5 LEW hepatocyte-specific signatures \n(based on set1 samples) over set2 samples UMAP')


### is varimax 15 strain specific
ggplot(set2.umap.df, aes( x=strain,y=var15head_UCell,fill=strain))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme(text = element_text(size=18))+
  ylab('Enrichment score of Lew signature in set2 map macrophages\n(based on set-1 Varimax-15)')

ggplot(set2.umap.df, aes( x=strain,y=var15tail_UCell,fill=strain))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = c("#CC79A7","#0072B2"))+ theme(text = element_text(size=22))

### is varimax 5 strain specific
ggplot(set2.umap.df, aes( x=strain,y=var5head_UCell,fill=strain))+geom_boxplot()+
  theme_classic()

ggplot(set2.umap.df, aes( x=strain,y=var5tail_UCell,fill=strain))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))


table(set2.umap.df$label[set2.umap.df$var5head_UCell>0.3])




set2.umap.df$cluster = as.character(set2_data$cluster)
#### checking if varimax-5 is valid in the Hep population of the set2 map
set2.umap.df.sub = set2.umap.df[set2.umap.df$cluster%in% c('8', '13', '11'),]
ggplot(set2.umap.df.sub, aes( x=strain,y=var15head_UCell,fill=strain))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme(text = element_text(size=13))+
  ggtitle('Enrichment score of Lewis signatures\nin set-2 map macrophage clusters\n(based on set-1 Varimax-15)')+
  ylab('Enrichment score of\nvarimax-15 Lewis signitures')+xlab('Strain')+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"


ggplot(set2.umap.df.sub, aes( x=strain,y=var15tail_UCell,fill=strain))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme(text = element_text(size=13))+
  ggtitle('Enrichment score of DA signatures\nin set-2 map macrophage clusters\n(based on set-1 Varimax-15)')+
  ylab('Enrichment score of\nvarimax-15 DA signitures')+xlab('Strain')+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"


#### checking if varimax-5 is valid in the Hep population of the set2 map
set2.umap.df.sub = set2.umap.df[str_sub(set2.umap.df$label, 1,3)=='Hep',]
ggplot(set2.umap.df.sub, aes( x=strain,y=var5head_UCell,fill=strain))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme(text = element_text(size=13))+
  ggtitle('Enrichment score of DA signatures\nin set-2 map hepatocyte clusters\n(based on set-1 Varimax-5)')+
  ylab('Enrichment score of\nvarimax-5 DA signitures')+xlab('Strain')+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"

ggplot(set2.umap.df.sub, aes( x=strain,y=var5tail_UCell,fill=strain))+geom_boxplot()+
  theme_classic()+scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme(text = element_text(size=13))+
  ggtitle('Enrichment score of Lewis signatures\nin set-2 map hepatocyte clusters\n(based on set-1 Varimax-5)')+
  ylab('Enrichment score of\nvarimax-5 Lewis signitures')+xlab('Strain')+
  stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"





###### checking consistency with the DE genes:

get_head_10genes(varimax_15_set1)[get_head_10genes(varimax_15_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes
get_tail_10genes(varimax_15_set1)[get_tail_10genes(varimax_15_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes

get_head_10genes(varimax_5_set1)[get_head_10genes(varimax_5_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes
get_tail_10genes(varimax_5_set1)[get_tail_10genes(varimax_5_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes

get_head_10genes(varimax_5_set1)[get_head_10genes(varimax_5_set1) %in% LEWvsDA_genes.df$Var1] ### should not give any results
get_tail_10genes(varimax_5_set1)[get_tail_10genes(varimax_5_set1) %in% DAvsLEW_genes.df$Var1] ### should not give any results




####### choosing the final genes for generation of dotplots ####### 
###### genes for Mac related strain differences - based on varimax-15
x=get_head_genes(varimax_15_set1)[get_head_genes(varimax_15_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes
y=get_tail_genes(varimax_15_set1)[get_tail_genes(varimax_15_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes
final_strain_markers = c(x,y) ### Lew first - DA second
final_strain_markers = final_strain_markers[!final_strain_markers %in% 'RT1-CE10'] 

###### genes for Hep related strain differences - based on varimax-5
x=get_head_genes(varimax_5_set1)[get_head_genes(varimax_5_set1) %in% DAvsLEW_genes.df$Var1] ### validated DA genes
y=get_tail_genes(varimax_5_set1)[get_tail_genes(varimax_5_set1) %in% LEWvsDA_genes.df$Var1] ### validated LEW genes
final_strain_markers = c(y,x) ### Lew first - DA second
final_strain_markers = final_strain_markers[!final_strain_markers %in% c('Apoc1', 'Scd')]


final_strain_markers =  c("Itgal","AC114233.2","Timp2","RT1-A1","Tmem176a" ,
                          "RT1-A2" ,"RGD1562667","Siglec5","RGD1559482","Ly6al")
#### only finding the strain genes based on varimax results
final_strain_markers = c(get_head_10genes(varimax_15_set1), 
                         get_tail_10genes(varimax_15_set1))

final_strain_markers = c(get_head_genes(varimax_15_set1), 
                         get_tail_genes(varimax_15_set1))


DoHeatmap((set2_data), features = final_strain_markers,slot='counts')

final_strain_markers_test = c(final_strain_markers, 'Icam1')


Idents(set1_data) = set1_data$strain
DotPlot(set1_data,features =final_strain_markers) + 
  RotatedAxis()+ xlab('Markers')+ylab('')+theme(axis.text.x = element_text(size=11))

table(set2_data$label)
table(set1_data$label)

### set2 subsets
set2_data_tmp = set2_data[,set2_data$cluster %in% as.character(c(17,10,11,8,12,7,13))] # all immune
set2_data_tmp = set2_data[,set2_data$cluster %in% as.character(c(11, 8, 13))] ## Mac
set2_data_tmp = set2_data[,set2_data$cluster %in% as.character(c(0, 1, 2, 3, 9, 14, 15))] ## Hep

## set1: Mac: 9, 10, 5
## set1: immune: 9, 10, 5, 13
set1_data_tmp = set1_data[,set1_data$cluster %in% as.character(c(9, 10, 5, 13))] # all immune
set1_data_tmp = set1_data[,set1_data$cluster %in% as.character(c(9, 10, 5))] ## Mac
Idents(set1_data_tmp) = set1_data$strain[set1_data$cluster %in% as.character(c(9, 10, 5))]

DotPlot(set1_data_tmp,features =final_strain_markers) + 
  RotatedAxis()+ xlab('Markers')+ylab('')+theme(axis.text.x = element_text(size=11))


Idents(set2_data_tmp) = set2_data_tmp$label

DotPlot(set2_data_tmp, split.by = "strain",
        features = final_strain_markers_test) + 
  RotatedAxis()+ xlab('Markers')+ylab('')+
  theme(axis.text.x = element_text(size=11))
#split.by = "groups"


############################ 
##### Frequency of lyz2 and Cd68 markers over CD45+ cells for each strain - check both set-1 and set2 maps
set1_data@meta.data$strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2))

data = set2_data
marker.df = data.frame(strain= data$strain,
                       Ptprc = ifelse(GetAssayData(data)[rownames(data) == 'Ptprc',] > 0, 'Ptprc+', 'Ptprc-'),
                       Cd68 = ifelse(GetAssayData(data)[rownames(data) == 'Cd68',] > 0,'Cd68+', 'Cd68-'),
                       Lyz2 = ifelse(GetAssayData(data)[rownames(data) == 'Lyz2',] > 0,'Lyz2+', 'Lyz2-'),
                       Marco = ifelse(GetAssayData(data)[rownames(data) == 'Marco',] > 0,'Macro+', 'Macro-'))


marker.df$Markers = paste0(marker.df$Ptprc, marker.df$Cd68, marker.df$Lyz2)
marker.df2 = data.frame(table(Strain=marker.df$strain, Markers=marker.df$Markers))
head(marker.df2)
ggplot(data=marker.df2, aes(x=Markers, y=Freq, fill=Strain)) +
  geom_bar(color='black', stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("pink1", "skyblue1"))+ylab('Number of cell in the Set-1 map')+
  theme_classic()+xlab('')+theme(axis.text.x = element_text(angle = 90, size=12))



#### Looking at Ptprc+ cells only 
marker.df3 = marker.df[marker.df$Ptprc=='Ptprc+',]

### evaluating double marker expression
marker.df3$Markers = paste0( marker.df3$Cd68, marker.df3$Lyz2)
marker.df4 = data.frame(table(Strain=marker.df3$strain, Markers=marker.df3$Markers))

### evaluating single marker expression
marker.df4 = data.frame(table(Strain=marker.df3$strain, Markers=marker.df3$Marco))

ggplot(data=marker.df4, aes(x=Markers, y=Freq, fill=Strain)) +
  geom_bar(color='black',stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  ylab('Number of Ptprc+ cell in the Set-1 map')+
  theme_classic()+xlab('')+
  theme(axis.text.x = element_text(angle = 90, size=12))


###################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


marker.df_bin = marker.df
marker.df_bin$sample = unlist(lapply(str_split(set2_data$umi, '_'), function(x) paste0(x[2:5], collapse = '_')))

marker.df_bin$Ptprc.Cd68 = paste0(marker.df_bin$Ptprc, marker.df_bin$Cd68)
marker.df_bin$Ptprc.Marco = paste0(marker.df_bin$Ptprc, marker.df_bin$Marco)
marker.df_bin$Ptprc.Lyz2 = paste0(marker.df_bin$Ptprc, marker.df_bin$Lyz2)

head(marker.df_bin)
marker.df_bin.l=data.frame(markers=c(marker.df_bin$Ptprc.Cd68, 
                                     marker.df_bin$Ptprc.Marco, 
                                     marker.df_bin$Ptprc.Lyz2),
                           strain = rep(marker.df_bin$strain,3),
                           sample = rep(marker.df_bin$sample,3))
table(marker.df_bin.l$markers)
marker.df_bin.l = marker.df_bin.l[marker.df_bin.l$markers %in% c('Ptprc+Cd68+', 'Ptprc+Lyz2+', 'Ptprc+Macro+'),]
head(marker.df_bin.l)

tab = table(marker.df_bin.l$markers,
            marker.df_bin.l$sample, 
            marker.df_bin.l$strain)
tab = data.frame(rbind(melt(tab[,,1]),melt(tab[,,2])))
colnames(tab) = c('markers','sample', 'freq')
tab$strain = unlist(lapply(str_split(tab$sample, '_'), '[[', 1))
tab = tab[tab$freq != 0,]

ggplot(data=tab, aes(x=markers, y=freq,fill=sample)) +
  geom_bar(color='black',stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  ylab('Cell count')+
  theme_classic()+xlab('')

tab2 = data_summary(tab[,-2], 'freq', c('markers','strain'))
ggplot(data=tab2, aes(x=markers, y=freq,fill=strain)) +
  geom_bar(color='black',stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  ylab('Cell count')+
  theme_classic()+xlab('')+
  geom_errorbar(aes(ymin=freq-sd, ymax=freq+sd), width=.2,
                position=position_dodge(.9)) 


#### set1 - significant relationship between strain and Cd68 in the Ptprc+ cells - p.val=1.215e-07
#### set1 - significant relationship between strain and Lyz2 in the Ptprc+ cells - p.val= 0.007681

#### set2 - significant relationship between strain and Cd68 in the Ptprc+ cells - p.val=1.067e-12
#### set2 - NO significant relationship between strain and Lyz2 in the Ptprc+ cells - p.val= 0.7605

conting_table = table(marker.df3$strain, marker.df3$Cd68) 
conting_table = table(marker.df3$strain, marker.df3$Lyz2)

fisher.test(conting_table)
chisq.test(conting_table)


###### check these differrences in cluster number-5 which includes the non-Inf Mac 
# do you think these differences might be related to varriations in number of cells ? sample imbalance
# if donw-sample, should we start from raw data? prenormalization might effect the results??


################ checking histology statistics
CD3_data = read.csv('Data/Histology/CD3.csv')
CD68_data = read.csv('Data/Histology/CD68.csv')
CD8_data = read.csv('Data/Histology/CD8.csv')

his_data = CD8_data
rownames(his_data) = paste0('layer_',his_data[,1])
his_data = his_data[,-1]
dim(his_data)
his_data_DA = his_data[,1:30]
his_data_LEW = his_data[,31:ncol(his_data)]

#### total differeces between the DA and LEW  
t.test(colSums(his_data_DA), colSums(his_data_LEW))
wilcox.test(colSums(his_data_DA), colSums(his_data_LEW))

#### CD3 general p-value: 0.03442 is significant
#### CD68 general p-value: 0.2374 is Not significant
#### CD8 general p-value: 0.3881 is Not significant

#### Layer-specific differeces between the DA and LEW ### alternative test
t.test.df = data.frame(layer=rownames(his_data), p.value=rep(0,10))
wilcoxon.test.df = data.frame(layer=rownames(his_data), p.value=rep(0,10))

for(i in 1:nrow(his_data_DA)){
  t.test.df[i,2] = t.test(his_data_DA[i,], his_data_LEW[i,])$p.value
  wilcoxon.test.df[i,2] = t.test(his_data_DA[i,], his_data_LEW[i,])$p.value
}
t.test.df$sig = t.test.df$p.value < 0.05/10
t.test.df$sig

wilcoxon.test.df$sig = wilcoxon.test.df$p.value < 0.05/10
wilcoxon.test.df$sig
#### CD3 layer-9 p-value: 0.002447611 is significant
#### CD68 not sig
#### CD8 not sig
