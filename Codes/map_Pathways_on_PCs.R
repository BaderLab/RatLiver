source('Codes/Functions.R')
Initialize()

#### Do NOT RUN THIS SECTION ####

# Basic function to convert mouse to rat gene names
convert_Mouse2rat_GeneList <- function(x){
  require("biomaRt")
  rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x , 
                   mart = mouse, 
                   attributesL = c("rgd_symbol"), 
                   martL = rat, uniqueRows=T)
  ratx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(ratx))
  return(ratx)
}

get_converted_geneSet <- function(a_geneSet){
  a_geneSet_genes <- geneIds(a_geneSet)
  a_geneSet_converted2rat <- convert_Mouse2rat_GeneList(a_geneSet_genes)
  return(a_geneSet_converted2rat)
}

gmtFile = '~/ChEA_2016.gmt'
geneSets <- getGmt(gmtFile)


### converting all the gene sets in the ChEA dataset
geneSets_converted2rat <- lapply(geneSets, get_converted_geneSet)
names(geneSets_converted2rat) <- names(geneSets)
length(geneSets_converted2rat)
geneSets_converted2rat <- readRDS('~/geneSets_converted2rat.rds')


## subsetting the genes-sets for your TF of interest
a_gene_set_name <- 'PPARG' # 'HNF4A' # GATA1
is_a_geneset <- sapply(str_split(rownames(data.frame(cbind(nGenes(geneSets)))),pattern = ' '),  
                       function(x) x[1]) == a_gene_set_name
sum(is_a_geneset)
matched_geneSets <- geneSets[is_a_geneset]

matched_geneSets_converted2rat <- lapply(matched_geneSets, get_converted_geneSet)
names(matched_geneSets_converted2rat) <- names(matched_geneSets)

matched_geneSets_converted2rat <- readRDS(paste0('~/', a_gene_set_name, '_geneSets_converted2rat.rds'))
PPARG_geneSets_converted2rat <- readRDS('~/PPARG_geneSets_converted2rat.rds')


geneSet_to_check <- lapply(matched_geneSets, geneIds)
names(geneSet_to_check) <- names(matched_geneSets)
df_geneset_elements <- data.frame(genes=unique(unlist(geneSet_to_check)))
upset_df <- data.frame(do.call(cbind, lapply(geneSet_to_check, 
                                             function(x) ifelse(df_geneset_elements$genes %in% x, 1, 0))))
rownames(upset_df) <- df_geneset_elements$genes
UpSetR::upset(upset_df, sets.bar.color = "#56B4E9",order.by = "freq")





#######################################
#### Calculating the enrichment of genesets using AUCell  ####

source('Codes/Functions.R')
Initialize()

geneSets_converted2rat <-  readRDS('~/geneSets_converted2rat.rds')
lapply(geneSets_converted2rat, length)
lapply(geneSets_converted2rat, function(x) paste0(sum(x%in% rownames(merged_samples)), ' , ' , length(x)))

### checking PPARG instead of all gene-sets
geneSets_converted2rat <- PPARG_geneSets_converted2rat
rat_geneSets <- sapply(1:length(geneSets_converted2rat), 
       function(i) GeneSet(geneSets_converted2rat[[i]], 
                           setName=names(geneSets_converted2rat)[i]), simplify =F)

names(rat_geneSets) = names(geneSets_converted2rat)


#### loading the original count matrix
load('~/XSpecies/for_scClustViz.RData')

merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')
exprMatrix <- as.matrix(GetAssayData(merged_samples))

rownames(exprMatrix)
getHead(exprMatrix)

cells_rankings <- AUCell_buildRankings(exprMatrix, 
                                       nCores=detectCores()-2, 
                                       plotStats=TRUE)

# cells_rankings <- readRDS('~/XSpecies/cells_rankings.rds')
saveRDS(cells_rankings, '~/XSpecies/cells_rankings_imputed.rds')
dim(cells_rankings)



AUCell_dir <- '~/XSpecies/Results/preproc_rats/merged/AUCell/'
###
cells_AUC <- sapply(1:length(rat_geneSets), function(i){ # 645 gene sets
  print('>>>>>>>>>>>>>')
  print(i)
  a_gene_set <- rat_geneSets[[i]]
  print(sum(geneIds(a_gene_set) %in% rownames(exprMatrix))/length(geneIds(a_gene_set)))
  a_gene_set_name <- gsub(' ', '_' ,names(rat_geneSets)[i])
  
  # the top 4 percent are being considered: can change with: aucMaxRank
  cells_AUC <- AUCell_calcAUC(a_gene_set, cells_rankings) 
  #saveRDS(cells_AUC, file=paste0(AUCell_dir, a_gene_set_name ,".rds"))
  return(cells_AUC)
}, simplify = F )




names(cells_AUC) <- names(rat_geneSets)
cells_AUC_df <- lapply(cells_AUC, function(x) data.frame(t(data.frame(getAUC(x)))))
cells_AUC_df <- do.call(cbind, cells_AUC_df)
cells_AUC_df$UMI=rownames(cells_AUC_df)

cells_AUC_df <- readRDS('~/cells_AUC_df.rds')
getHead(cells_AUC_df)
dim(cells_AUC_df)

### ???
PPARG_mac_AUC_obj <- readRDS('~/XSpecies/Results/preproc_rats/merged/AUCell/PPARG_20176806_ChIP-Seq_MACROPHAGES_Mouse.rds')
###

rot_data <- readRDS('Results/preproc_rats/merged/rot_data.rds')
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
scores$UMI <- rownames(scores)

# [1] "PPARG.19300518.ChIP.PET.3T3.L1.Mouse"        
# [2] "PPARG.20176806.ChIP.Seq.THIOMACROPHAGE.Mouse"
# [3] "PPARG.23326641.ChIP.Seq.C3H10T1.2.Mouse"     
# [4] "PPARG.20176806.ChIP.Seq.3T3.L1.Mouse"        
# [5] "PPARG.20176806.ChIP.Seq.MACROPHAGES.Mouse"   
# [6] "PPARG.20887899.ChIP.Seq.3T3.L1.Mouse" 

AUcell_PC_matrix <- merge(cells_AUC_df, scores, 'UMI', 'UMI', all.y=T, sort=F)
colnames(AUcell_PC_matrix)[2:7]

pdf('Results/preproc_rats/merged/AUCell/pc_results.pdf')
ggplot(AUcell_PC_matrix, 
       aes(x=PC_1, y=PC_5, color=PPARG.20887899.ChIP.Seq.3T3.L1.Mouse))+
  geom_point()+theme_bw()+scale_color_viridis(direction = -1)
dev.off()



saveRDS(merged_samples, 'Results/preproc_rats/merged/merged_rat_samples_2.rds')




