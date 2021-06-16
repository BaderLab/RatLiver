### annotation based on Tallulah markers - rds files generated in the clean_marker_set.R

Total_markers_converted_df <- readRDS('~/XSpecies/Data/Total_markers_converted_df.rds') ## all the markers merged together
immune_markers <- readRDS('~/XSpecies/Data/McParland_markers/liver_markers_tallulah/liver_immune_markers_mapped.rds') ## immune cell markers
general_markers_df <- readRDS('~/XSpecies/Data/McParland_markers/liver_markers_tallulah/liver_general_markers.rds') ## liver general markers

head(immune)


to_include <- function(items) which(items != '' & !is.na(items))


### getting the cleaned list of Tallulah markers
source('~/XSpecies/Codes/convert_human_to_ortholog_functions.R')


########################################################################
#### import input markers file: CellType markers - Liver ####
markers_liver <- read.csv('Markers/Celltype markers - Liver.csv')

### adding to human attributes
markers_liver_human <- markers_liver[markers_liver$Species=='Human',]
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

#### convert the human ensemble ids to rat orthologs #### 
markers_liver_human_mapped = .getMapped_hs2model_df(ensembl, 
                                                      candidateGenes=markers_liver_human$Gene, 
                                                      model_animal_name = "rnorvegicus")
markers_liver_human_mapped <- merge(markers_liver_human, markers_liver_human_mapped, by.x='Gene', by.y='symbol') 
markers_liver_human_mapped <- markers_liver_human_mapped[to_include(markers_liver_human_mapped$rnorvegicus_homolog_associated_gene_name),]
markers_liver_human_mapped.df <- data.frame(gene=markers_liver_human_mapped$rnorvegicus_homolog_associated_gene_name,
                                            cell_type=markers_liver_human_mapped$Celltype, 
                                            Species='Human')
head(markers_liver_human_mapped.df)


#### adding to mouse attributes ####
markers_liver_mouse <- markers_liver[markers_liver$Species=='Mouse',]
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('mmusculus_gene_ensembl',mart=ensembl)
#### convert the mouse ensemble ids to rat orthologs
markers_liver_mouse_mapped = .getMapped_mouse2model_df(markers_liver_mouse$Gene)
markers_liver_mouse_mapped <- merge(markers_liver_mouse, markers_liver_mouse_mapped, by.x='Gene', by.y='MGI.symbol') 
markers_liver_mouse.df = data.frame(gene=markers_liver_mouse_mapped$RGD.symbol,
                                    cell_type=markers_liver_mouse_mapped$Celltype,
                                    Species='Mouse') 
head(markers_liver_mouse.df)



### adding to rat attributes ### 
markers_liver_rat <- markers_liver[markers_liver$Species=='Rat',]
markers_liver_rat_df <- get_rat_ensembl_ids(markers_liver_rat$Gene)
markers_liver_rat_df <- merge(markers_liver_rat, markers_liver_rat_df, by.x='Gene', by.y='rgd_symbol', all.x=T)
markers_liver_rat.df <- data.frame(gene=markers_liver_rat_df$Gene,
                                   cell_type=markers_liver_rat_df$Celltype, 
                                   Species='Rat')  

#### merging all the 3 species markers ####
markers_liver_mapped <- rbind(markers_liver_human_mapped.df,markers_liver_mouse.df,markers_liver_rat.df)
markers_liver_mapped <- markers_liver_mapped[!duplicated(markers_liver_mapped), ]
markers_liver_mapped.list <- split(markers_liver_mapped, markers_liver_mapped$cell_type)


########################################################################



########################################################################
#### import input markers file: CellType markers - Immune ####
markers_immune <- read.csv('Markers/Celltype markers - Immune.csv', header = F)
head(markers_immune)
table(markers_immune$V3)
### adding to human attributes
markers_immune_human <- markers_immune[markers_immune$V3=='Human',]
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

#### convert the human ensemble ids to rat orthologs #### 
markers_immune_human_mapped = .getMapped_hs2model_df(ensembl, 
                                                    candidateGenes=markers_immune_human$V1, 
                                                    model_animal_name = "rnorvegicus")
markers_immune_human_mapped <- merge(markers_immune_human, markers_immune_human_mapped, by.x='V1', by.y='symbol') 
markers_immune_human_mapped <- markers_immune_human_mapped[to_include(markers_immune_human_mapped$rnorvegicus_homolog_associated_gene_name),]
markers_immune_human.df <- data.frame(gene=markers_immune_human_mapped$rnorvegicus_homolog_associated_gene_name,
                                            cell_type=markers_immune_human_mapped$V2, 
                                            Species='Human')
head(markers_immune_human.df)

#### adding to mouse attributes ####
markers_immune_mouse <- markers_immune[markers_immune$V3=='Mouse',]
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('mmusculus_gene_ensembl',mart=ensembl)
#### convert the mouse ensemble ids to rat orthologs
markers_immune_mouse_mapped = .getMapped_mouse2model_df(markers_immune_mouse$V1)
markers_immune_mouse_mapped <- merge(markers_immune_mouse, markers_immune_mouse_mapped, by.x='V1', by.y='MGI.symbol') 
markers_immune_mouse.df = data.frame(gene=markers_immune_mouse_mapped$RGD.symbol,
                                    cell_type=markers_immune_mouse_mapped$V2,
                                    Species='Mouse') 
head(markers_immune_mouse.df)

markers_immune_mapped <- rbind(markers_immune_human.df, markers_immune_mouse.df)
markers_immune_mapped <- markers_immune_mapped[!duplicated(markers_immune_mapped), ]

head(markers_immune_mapped)
#saveRDS(markers_immune_mapped, 'Markers/markers_immune_mapped.rds')

########################################################################

#### import input markers file: CellType markers - Liver ####
markers_subclustering <- read.csv('Markers/Celltype markers - SubClusteringGeneLists.csv')






########################################################################
###### checking the presence of well-known cell-type markers in the top DE genes

markers_liver_mapped.list <- split(markers_liver_mapped, markers_liver_mapped$cell_type)
markers_immune_mapped.list <- split(markers_immune_mapped, markers_immune_mapped$cell_type)

lapply(markers_liver_mapped.list, head)
lapply(markers_immune_mapped.list, head)

Cluster_markers <- readRDS('Results/Cluster_markers.rds')

top_gene_number = 30

for(i in 1:length(Cluster_markers)){
  print('------------------------------------------------------------')
  print(names(Cluster_markers)[i])
  print('------------------')
  a_cluster_DE_genes <- rownames(Cluster_markers[[i]])[1:top_gene_number]
  
  for(j in 1:length(markers_immune_mapped.list)){
    DE_genes_included <- a_cluster_DE_genes[a_cluster_DE_genes %in% markers_immune_mapped.list[[j]]$gene]
    if(length(DE_genes_included)>0) {
      print(names(markers_immune_mapped.list)[j])
      print(DE_genes_included)}
  }
}

#############################################################################
############### generating the ranked files for the pathway analysis ########

# smallest_number in R: 
smallest_number <- .Machine$double.xmin

for(i in 1:length(Cluster_markers)){
  a_cluster_name <- names(Cluster_markers)[i]
  a_cluster_markers <- Cluster_markers[[i]]
  a_cluster_markers.df <- data.frame(genes=rownames(a_cluster_markers),
                                     rank=-log10(a_cluster_markers$p_val_adj+smallest_number)*sign(a_cluster_markers$avg_log2FC))
  
  a_cluster_markers.df.sorted <- a_cluster_markers.df[order(a_cluster_markers.df$rank, decreasing = T),]
  write.table(a_cluster_markers.df.sorted, 
              paste0('Results/ranked_files/',a_cluster_name, '.rnk'),
              quote = F, col.names = F, row.names = F, sep = '\t')
}






