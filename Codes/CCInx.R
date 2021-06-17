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
