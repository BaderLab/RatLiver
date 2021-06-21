source('Codes/Functions.R')
Initialize()
library(CCInx)

split_inx.edge <- function(inx, aStrain){
  inx_edge <- inx$edges[order(inx$edges$edgeWeight, decreasing=T),]
  cellTypeA = sapply(str_split(inx_edge$nodeA, '_'), '[[', 2)
  cellTypeB = sapply(str_split(inx_edge$nodeB, '_'), '[[', 2)
  colnames(inx_edge)[3] = paste0(colnames(inx_edge)[3], '_', aStrain)
  inx_edge$cellTypes = paste0(cellTypeA, '_', cellTypeB)
  inx_edge$nodeAB = paste(inx_edge$nodeA, inx_edge$nodeB)
  inx_edge.split <- split.data.frame(inx_edge, inx_edge$cellTypes)
  return(inx_edge.split)
}

############ loading the data ############
#old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
load(new_data_scCLustViz_object)


#### strain comparison ####
strains <- sapply(str_split(colnames(your_scRNAseq_data_object),'_'), '[[', 2)
aStrain = 'LEW' #'LEW'#"DA"

#### generating the meta dataframe ####
annot_df <- readRDS('cluster_labels/new_labels.rds')
annot_df$label <- paste0(annot_df$label,'(',annot_df$cluster,')')
annot_df <- annot_df[,c('cluster', 'label')]
head(annot_df)

meta = data.frame(labels=sCVdata_list$res.0.6@Clusters, 
                  umi=names(sCVdata_list$res.0.6@Clusters),
                  strain=sapply(str_split(names(sCVdata_list$res.0.6@Clusters),'_'), '[[', 2))
head(meta)

merged_meta <- merge(meta, annot_df, by.x='labels', by.y='cluster', all.x=T,sort=F)
merged_meta <- merged_meta[match(meta$umi,merged_meta$umi),]
head(merged_meta)

### finding the list of genes to be removed for lack of orthology
rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
mouse_ortho_genes <- rat_to_mouse_genes$mmusculus_homolog_associated_gene_name[match(rownames(your_scRNAseq_data_object),
                                                                                     rat_to_mouse_genes$symbol)]
to_be_removed <- is.na(mouse_ortho_genes) | mouse_ortho_genes==''
sum(to_be_removed) ## 1522 genes removed


#### defining the input sample ####
data.input <- GetAssayData(your_scRNAseq_data_object)
rownames(data.input) <- make.unique(mouse_ortho_genes)
data.input <- data.input[!to_be_removed,]

#### strain comparison ####
data.input <- data.input[,strains==aStrain]
meta <- meta[meta$strain==aStrain,]
merged_meta <-  merged_meta[merged_meta$strain==aStrain,]

data.input.seur <- CreateSeuratObject(counts = data.input, 
                                      row.names = rownames(data.input), 
                                      assay = "logcounts") 

gsl <- BuildGeneStatList(inD=data.input.seur,
                         cl=merged_meta$label,
                         assayType="logcounts")
lapply(gsl[1:3],head)

inx <- BuildCCInx(GeneStatList=gsl,
                  Species="mmusculus")

head(inx$edges)
head(inx$nodes)

saveRDS(inx, paste0('Objects/CCinx_', aStrain, '.rds'))
inx <- readRDS(paste0('Objects/CCinx_', aStrain, '.rds'))

inx_DA <- readRDS('Objects/CCinx_DA.rds')
inx_LEW <- readRDS('Objects/CCinx_LEW.rds')

cell_types_interest <-  c("Inf Mac(11)", "LSEC(3)", "LSEC(4)", "LSEC(6)", "Mac(13)",
                          "Mac(8)","Non-Inf Mac(17)", "Plasma(12)","Stellate(16)", "T cell(7)")




inx_edge.split.DA <- split_inx.edge(inx_DA, 'DA')
inx_edge.split.LEW <-split_inx.edge(inx_LEW, 'LEW') 

names(inx_edge.split.DA) == names(inx_edge.split.LEW)

inx_edge.split.merged = sapply(1:length(inx_edge.split.DA), 
       function(i){
         merge.df = merge(inx_edge.split.DA[[i]], inx_edge.split.LEW[[i]], by.x='nodeAB', by.y='nodeAB')
         merge.df = merge.df[,!colnames(merge.df) %in% c('nodeA.y', 'nodeB.y', 'cellTypes.y') ]
         merge.df$WeightDiff = (merge.df$edgeWeight_LEW) - (merge.df$edgeWeight_DA)
         merge.df <- merge.df[order(merge.df$WeightDiff),]
         return(merge.df)
         }, simplify = F)

names(inx_edge.split.merged) = names(inx_edge.split.DA)
lapply(inx_edge.split.merged, head)
hist(inx_edge.split.merged[[1]]$WeightDiff)

cellTypeA = 'LSEC(3)'
cellTypeB = 'Mac(13)'
LR_A = 'Fn1' #'Fn1' 'Lamc1' Vtn
LR_B = 'Itgb1'

ABintX <- inx_edge.split.merged[[paste0(cellTypeA, '_', cellTypeB)]]
hist(ABintX$WeightDiff)
ABintX[ABintX$nodeA.x == paste0(LR_A, '_', cellTypeA) & ABintX$nodeB.x == paste0(LR_B, '_', cellTypeB),] 



inx_edge.split.DA <- inx_edge.split
head(inx_LEW$edges)

pdf('Plots/set2map_CCIInx_strainComparison.pdf',height = 20, width=25)
for (i in 1:length(cell_types_interest)){
  for (j in i:length(cell_types_interest)){
    par(mfrow=c(1,2))
    cellTypeA = cell_types_interest[i]
    cellTypeB = cell_types_interest[j]
      
    PlotCCInx(INX=inx_DA,
              cellTypeA=cellTypeA,
              cellTypeB=cellTypeB,
              proteinTypeA="Ligand",
              proteinTypeB="Receptor",
              GeneMagnitudeThreshold=1)
    
    PlotCCInx(INX=inx_LEW,
              cellTypeA=cellTypeA,
              cellTypeB=cellTypeB,
              proteinTypeA="Ligand",
              proteinTypeB="Receptor",
              GeneMagnitudeThreshold=1)
    
  }
}

dev.off()


######################################################
####### evaluating ccx based on varimax factors - set-1 map 
####### we would like to evaluate whether Mac<->LSECs interactions are more enriched in one strain
##### varimax-4 represents cluster-3, annotated as central venous LSECs+portal endothelial cells 
##### -(varimax-7) represent cluster 7, 14, annotated as activated LSECs/Hep stellate cells
##### varimax-15 positive scores represent LEW-enriched genes and neg values represent DA-enriched


### importing the data and converting rat gene-names to mouse orthologs
varimax_res <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
mouse_ortho_genes <- rat_to_mouse_genes$mmusculus_homolog_associated_gene_name[match(rownames(varimax_res$rotLoadings),
                                                                                     rat_to_mouse_genes$symbol)]
to_be_removed <- is.na(mouse_ortho_genes) | mouse_ortho_genes==''
sum(to_be_removed) ## 1522 genes removed

varimax_load <- data.frame(sapply(1:ncol(varimax_res$rotLoadings),
                       function(i) varimax_res$rotLoadings[,i]))
colnames(varimax_load) <- paste0('var', colnames(varimax_res$rotLoadings))

varimax_load <- varimax_load[!to_be_removed,]
row.names(varimax_load) <- make.unique(mouse_ortho_genes[!to_be_removed])
#row.names(varimax_load) <- str_replace(rownames(varimax_load), '_', '-')

#### adding the reversed directions of varimax-pc 7 and 15
varimax_load$varPC_15Rev <- -(varimax_load$varPC_15)
varimax_load$varPC_7Rev <- -(varimax_load$varPC_7)
varimax_load$varPC_5Rev <- -(varimax_load$varPC_5)

#### varimax factors to look into

#### Macrophage_PCs
group1 = data.frame(rbind(c('varPC_6', 'cluster-9', 'Inf-Mac'), 
              c('varPC_1', 'cluster-10', 'nonInf-Mac'), 
              c('varPC_21', 'cluster-5', 'nonInf-Mac'),
              ## strain-specific Mac PCs
              c('varPC_15', 'cluster-5,14,10', 'LEW-nonInf-Mac'), 
              c('varPC_15Rev', 'cluster-5,14,10', 'DA-nonInf-Mac')))




#### T/NK cells
group2 = data.frame(rbind(c('varPC_10', 'cluster-13', 'T-NK-cells'),
              #### LSECs PCs
              c('varPC_4', 'cluster-3', 'cv-LSECs'), 
              c('varPC_7Rev', 'cluster-7,14', 'LSECs')))


colnames(group1) = c('varPC', 'cluster', 'cell-type')
colnames(group2) = c('varPC', 'cluster', 'cell-type')

group1$group = 'group1'
group2$group = 'group2'
group = rbind(group1, group2)
group$names = paste0(group$varPC,sep = ' ; ', group$cluster,sep = ' ; ', group$`cell-type` )

### all cells in a strain
group3 = list(c('varPC-5', 'all-clusters', 'LEW-all'),
              c('varPC-5Rev', 'all-clusters', 'DA-all'))


factors_list <- c(group1$varPC, group2$varPC)

#### generating a list of dataframes similar to DE list from the varimax factors 
varimax_load.list <- sapply(1:ncol(varimax_load),
                       function(i) data.frame(logFC=varimax_load[,i], 
                                              padj=0.01, 
                                              row.names = rownames(varimax_load)), 
                       simplify = F)

names(varimax_load.list) <- colnames(varimax_load)

varimax_load.list <- varimax_load.list[names(varimax_load.list) %in% factors_list]
names(varimax_load.list) <- sapply(names(varimax_load.list), function(a_name) group$names[group$varPC == a_name])
names(varimax_load.list) <- str_replace(names(varimax_load.list), '_', '')
group$names <- str_replace(group$names, '_', '')

inx <- BuildCCInx(GeneStatList=varimax_load.list,
                  GeneMagnitude="logFC",
                  Species="mmusculus")

par(mfrow=c(3,3))
for(i in 1:(length(varimax_load.list))){
  load.mat <- varimax_load.list[[i]]
  hist(load.mat$logFC,breaks = seq(-0.3, 0.3, 0.02), main=names(varimax_load.list)[i], 
       col=brewer.pal(n = 8, name = "Set1")[i], #RdBu
       xlab='varimax-score')
}

sapply(1:length(varimax_load.list), 
       function(i){
         load.mat <- varimax_load.list[[i]]
         sum(load.mat>0.07)
         })



is.group1 = group$group == 'group1'
is.group2 = group$group == 'group2'

pdf('Plots/set1map_CCInx_varimaxBased_2.pdf',height = 10, width=25)

for(i in 1:sum(is.group1)){
  for(j in 1:sum(is.group2)){
    
    par(mfrow=c(1,2))
    PlotCCInx(INX=inx,
              cellTypeA=group$names[is.group1][i],
              cellTypeB=group$names[is.group2][j],
              proteinTypeA="Ligand",
              proteinTypeB="Receptor",
              TopEdges=60)
    
    PlotCCInx(INX=inx,
              cellTypeA=group$names[is.group2][j],
              cellTypeB=group$names[is.group1][i],
              proteinTypeA="Ligand",
              proteinTypeB="Receptor",
              TopEdges=60)
  }
}
dev.off()




pdf('Plots/set2map_CCIInx_strainComparison.pdf',height = 20, width=25)
for (i in 1:length(cell_types_interest)){
  for (j in i:length(cell_types_interest)){
    par(mfrow=c(1,2))
    cellTypeA = cell_types_interest[i]
    cellTypeB = cell_types_interest[j]
    
    PlotCCInx(INX=inx_DA,
              cellTypeA=cellTypeA,
              cellTypeB=cellTypeB,
              proteinTypeA="Ligand",
              proteinTypeB="Receptor",
              GeneMagnitudeThreshold=1)
    
    PlotCCInx(INX=inx_LEW,
              cellTypeA=cellTypeA,
              cellTypeB=cellTypeB,
              proteinTypeA="Ligand",
              proteinTypeB="Receptor",
              GeneMagnitudeThreshold=1)
    
  }
}

dev.off()
