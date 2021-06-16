# devtools::install_github("sqjin/CellChat")
library(caret)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#### loading the data
#old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"
load(new_data_scCLustViz_object)

### finding the list of genes to be removed for lack of orthology
rat_to_mouse_genes <- readRDS('~/XSpecies/rat_to_mouse_genes.rds')
mouse_ortho_genes <- rat_to_mouse_genes$mmusculus_homolog_associated_gene_name[match(rownames(your_scRNAseq_data_object),
                                                                                     rat_to_mouse_genes$symbol)]
to_be_removed <- !(is.na(mouse_ortho_genes) | mouse_ortho_genes=='')
sum(!to_be_removed) ## 1522 genes removed

#### loading the annotations
#annot_df <- readRDS('cluster_labels/old_labels.rds')
annot_df <- readRDS('cluster_labels/new_labels.rds')
annot_df$label <- paste0(annot_df$label,'(',annot_df$cluster,')')
annot_df <- annot_df[,c('cluster', 'label')]

#### strain comparison ####
strains <- sapply(str_split(colnames(your_scRNAseq_data_object),'_'), '[[', 2)
aStrain = 'DA' #'LEW'#"DA"

#### generating the meta dataframe ####
meta = data.frame(labels=sCVdata_list$res.0.6@Clusters, 
                  umi=names(sCVdata_list$res.0.6@Clusters),
                  strain=sapply(str_split(names(sCVdata_list$res.0.6@Clusters),'_'), '[[', 2))

#### down-sampling #### 
meta_split <- split.data.frame(meta, meta$labels)
meta_split_ds <- lapply(meta_split, function(df) {downSample(df, as.factor(df$strain), 
                                                             list = FALSE, 
                                                             yname = "Class")})
meta_ds <- do.call(rbind, meta_split_ds) ## down-sampled meta dataframe
meta <- meta_ds


#### defining the input sample ####
data.input <- GetAssayData(your_scRNAseq_data_object)
rownames(data.input) <- make.unique(mouse_ortho_genes)

#### down-sampling the meta-data frame and the data.input ####
data.input <- data.input[to_be_removed,(colnames(data.input) %in% meta$umi) & strains==aStrain]
#### strain comparison ####
meta <- meta[meta$strain==aStrain,]

##### checking the dimensions of input data and meta dataframe
dim(data.input)
dim(meta)

merged_meta <- merge(meta, annot_df, by.x='labels', by.y='cluster', all.x=T,sort=F)
merged_meta <- merged_meta[match(meta$umi,merged_meta$umi),]
head(merged_meta)

unique(meta$labels) # check the cell labels
cellchat <- createCellChat(object = data.input, meta = merged_meta, group.by = "label")

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)


# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multiprocess", workers = detectCores()-2) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
unique.data.frame(df.net[,c('ligand','receptor')])
df.net2 <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat, paste0('Objects/cellchat_', aStrain, '.rds'))



cellchat <- readRDS(paste0('Objects/cellchat_', aStrain, '.rds'))


par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat@DB

### visualizing the predicted CCC network ####
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 13:18) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#### evaluating the enriched pathways####
cellchat@netP$pathways
length(cellchat@netP$pathways)

par(mfrow = c(2,3), xpd=TRUE)
for (i in 13:15) {
  pathways.show <- cellchat@netP$pathways[i]
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
}

#### visualizing one pathway 
pathways.show <- c("GALECTIN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")



#################################################
########## Comparative analysis of CCC in the LEW and DA rats ####### 

cellchat.DA <- readRDS('Objects/cellchat_DA.rds')
cellchat.LEW <- readRDS('Objects/cellchat_LEW.rds')
object.list <- list(DA = cellchat.DA, LEW = cellchat.LEW)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

## is it required to apply liftCellChat when the cells from both populations are present in all the clusters?
##
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1+gg2


# red (or blue) colored edges represent increased (or decreased) 
# signaling in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#  red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


#### Comparing the outgoing and incoming interaction strength in 2D space ####
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


#### comparing the counts of interaction in each cell type ####
#### directly show the number of interactions or interaction strength between
#### any two cell populations in each dataset

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#### comparing the signaling pathways between the two groups ####
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2



##### understanding the output data struture
cellchat.DA <- readRDS('Objects/cellchat_DA.rds')

dim(cellchat.DA@net$pval)
dim(cellchat.DA@net$count)
dim(cellchat.DA@net$prob)

names(cellchat.DA@net$pval[,1,1])
dim(cellchat.DA@net$pval[1,,])
cellchat.DA@net$pval[1,1,]

############################
#### Implementing the cellchat.DA@net$count matrix 
#### also, finding the list of significant LR for each cell group pairs 

dim1_names <- names(cellchat.DA@net$pval[,1,1])
dim2_names <- names(cellchat.DA@net$pval[1,,1])
dim3_names <- names(cellchat.DA@net$pval[1,1,])

ccx_mat_pval  = data.frame(matrix(data = rep(-1, length(dim1_names)^2), 
                                  nrow = length(dim1_names), 
                                  ncol = length(dim2_names)))
rownames(ccx_mat_pval) = dim1_names
colnames(ccx_mat_pval) = dim2_names
LR_list = list()
LR_list_prob = list()

for(dim1 in dim1_names){
  for(dim2 in dim2_names){
    cell_cell_p_val = cellchat.DA@net$pval[dim1, dim2, ]
    
    ### filing the count matrix (#sig interaction between cell pairs)
    ccx_mat_pval[dim1, dim2] = sum(cell_cell_p_val<0.05)
    
    ### filing the LR list (name of sig interaction between cell pairs)
    LR_list[[dim1]][dim2] = list(names(cell_cell_p_val[cell_cell_p_val<0.05]))
    
    ### filing the LR/prob list (name of sig interaction between cell pairs and their probabilities)
    if(length(LR_list[[dim1]][[dim2]])>0)
      LR_list_prob[[dim1]][[dim2]] = cellchat.DA@net$prob[dim1, dim2,][LR_list[[dim1]][[dim2]]]
   
  }
}


pheatmap(cellchat.DA@net$count)
pheatmap(ccx_mat_pval)
############################
LR_list$`LSEC(3)`
LR_list_prob$`LSEC(3)`

names(cell_cell_p_val[cell_cell_p_val<0.05])


dim1 = dim1_names[1]
dim2 = dim2_names[1]
cell_cell_p_val = cellchat.DA@net$pval[dim1, dim2, ]
cell_cell_p_val[cell_cell_p_val<0.05]


cellchat.DA@var.features$features.info
head(cellchat.DA@DB$interaction)
head(cellchat.DA@DB$geneInfo)

head(cellchat.DA@LR$LRsig)

head(cellchat.DA@net$prob)
head(cellchat.DA@netP$prob)
head(cellchat.DA@net$pval)

?computeCommunProb



