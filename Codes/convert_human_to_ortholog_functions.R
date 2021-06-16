## Attension:
## you need to set the animal name here before loading it in another script

source('Codes/Functions.R')
Initialize()

model_animal_name = 'rnorvegicus' # 'mmusculus'

listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

.getMapped_hs2model_df <- function(ensembl=ensembl, candidateGenes, model_animal_name){
  
  wanted_attributes <- c(paste0(model_animal_name, '_homolog_ensembl_gene'), 
                         paste0(model_animal_name, '_homolog_associated_gene_name'), 
                         paste0(model_animal_name, '_homolog_orthology_type'),
                         paste0(model_animal_name, '_homolog_perc_id'), 
                         paste0(model_animal_name, '_homolog_perc_id_r1'), 
                         #paste0(model_animal_name, '_homolog_dn'), 
                         #paste0(model_animal_name, '_homolog_ds') , 
                         paste0(model_animal_name, '_homolog_orthology_confidence'))
  
  
  mappedGenesToOrthologs <- getBM(filters="hgnc_symbol", 
                                  attributes= c("ensembl_gene_id", wanted_attributes),
                                  values=candidateGenes, mart= ensembl)
  ## cleaning the resulting data frame
  hgnc_symbol_to_ensembl <- data.frame(getBM(filters="hgnc_symbol", 
                                             attributes= c("ensembl_gene_id", 'hgnc_symbol'),
                                             values=candidateGenes, mart= ensembl))
  print(dim(mappedGenesToOrthologs))
  print(length(candidateGenes))
  print(dim(hgnc_symbol_to_ensembl))
  
  candidateGenes.df <- data.frame(symbol=candidateGenes)
  candidateGenes.df <- merge(candidateGenes.df, hgnc_symbol_to_ensembl, by.x='symbol', by.y='hgnc_symbol',all.x=T)
  mappedGenesToOrthologs <- merge(candidateGenes.df, mappedGenesToOrthologs, by.x='ensembl_gene_id', by.y='ensembl_gene_id', all.x=T )
  
  return(mappedGenesToOrthologs)
}
# usage:
# .getMapped_hs2model_df(ensembl, candidateGenes, model_animal_name)





#  'rnorvegicus'  'mmusculus'
model_animal_name = 'mmusculus'

listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('rnorvegicus_gene_ensembl',mart=ensembl)

.getMapped_rat2model_df <- function(ensembl, candidateGenes, model_animal_name){
  
  wanted_attributes <- c(paste0(model_animal_name, '_homolog_ensembl_gene'), 
                         paste0(model_animal_name, '_homolog_associated_gene_name'), 
                         paste0(model_animal_name, '_homolog_orthology_type'),
                         paste0(model_animal_name, '_homolog_perc_id'), 
                         paste0(model_animal_name, '_homolog_perc_id_r1'), 
                         paste0(model_animal_name, '_homolog_orthology_confidence'))
  
  
  mappedGenesToOrthologs <- getBM(filters="rgd_symbol", 
                                  attributes= c(wanted_attributes, "ensembl_gene_id"),
                                  values=candidateGenes, mart= ensembl)
  ## cleaning the resulting data frame
  rgd_symbol_to_ensembl <- data.frame(getBM(filters="rgd_symbol", 
                                             attributes= c("ensembl_gene_id", 'rgd_symbol'),
                                             values=candidateGenes, mart= ensembl))
  print(dim(mappedGenesToOrthologs))
  print(length(candidateGenes))
  print(dim(rgd_symbol_to_ensembl))
  
  candidateGenes.df <- data.frame(symbol=candidateGenes)
  candidateGenes.df <- merge(candidateGenes.df, rgd_symbol_to_ensembl, by.x='symbol', by.y='rgd_symbol',all.x=T,sort=F)
  mappedGenesToOrthologs <- merge(candidateGenes.df, mappedGenesToOrthologs, by.x='ensembl_gene_id', by.y='ensembl_gene_id', all.x=T, sort=F)
  
  return(mappedGenesToOrthologs)
}



  
# Basic function to convert mouse to rat gene names
.getMapped_mouse2model_df <- function(x){
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
  return(genesV2)
}






#### #### #### #### converting the list of genes to their accourding orthologs
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)
ensembl = useDataset('rnorvegicus_gene_ensembl',mart=ensembl)


## checking the appropiate filter and attribute to use
listFilters(ensembl)[grep(listFilters(ensembl)[,1], pattern = 'symbol'),]
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'rnorvegicus'),] 
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'mmusculus'),] 
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'symbol'),] 
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'hgnc'),] 




####  Convert human gene list to ensembl id
get_human_ensembl_ids <- function(human_gene_symbol_list){
  listMarts()
  ensembl <- useMart("ensembl")
  datasets <- listDatasets(ensembl)
  ensembl = useDataset('hsapiens_gene_ensembl',
                       mart=ensembl)
  general_markers_human_ensembl_df <- getBM(filters="hgnc_symbol", 
                                            attributes= c('hgnc_symbol',"ensembl_gene_id"),
                                            values=as.character(human_gene_symbol_list), mart= ensembl)
  return(general_markers_human_ensembl_df)
}

get_rat_ensembl_ids <- function(rat_gene_symbol_list){
  listMarts()
  ensembl <- useMart("ensembl")
  datasets <- listDatasets(ensembl)
  ensembl = useDataset('rnorvegicus_gene_ensembl',
                       mart=ensembl)
  general_markers_rat_ensembl_df <- getBM(filters="rgd_symbol", 
                                            attributes= c('rgd_symbol',"ensembl_gene_id"),
                                            values=as.character(rat_gene_symbol_list), mart= ensembl)
  return(general_markers_rat_ensembl_df)
}



convertRatGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("rgd_symbol"), 
                   filters = "rgd_symbol", 
                   values = x , mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(genesV2) # return(humanx)
}

# convert rat to human
convert_Rat2human_GeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("rgd_symbol"), 
                   filters = "rgd_symbol", 
                   values = x , mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(genesV2) # return(humanx)
}


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


# Basic function to convert human to mouse gene names
convert_Human2mouse_GeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


