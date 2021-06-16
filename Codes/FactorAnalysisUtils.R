
elbow_plot <- function(standardDev, title=''){
  ## Extract the standard deviation value for each Principal component and plot them
  ## INPUTS:
  ##  seurat_obj: seurat object
  ## RETURNS:
  ##  elbow (scree) plot for the dimension reduction
  
  ndims = 50
  percVar = (standardDev^2 * 100)/sum(standardDev^2)
  plot <- ggplot(data = data.frame(dims = 1:ndims, percVar = percVar[1:ndims])) +
    geom_point(mapping = aes_string(x = 'dims', y = 'percVar')) +
    labs( x= 'PCs', y = 'Percentage of explained variance')+ggtitle(title)
  theme_classic()
  print(plot)
}



convert_features <- function(sample){
  sample_data <- GetAssayData(sample, assay = 'RNA')
  rownames(sample_data) <- mapper$V2[match(rownames(sample_data), mapper$V1)]
  return(CreateSeuratObject(sample_data))
}


get_varimax_rotated <- function(gene_exp_matrix, loading_matrix){
  
  ## gene_exp_matrix: gene expression matrix. rows named as genes and columns named as UMIs.
  ##                  attention: Data SHOULD be SCALED.
  ##                  Potentially gained from Seurat GetAssayData() function
  ## loading_matrix: the PCA loading matrix with rows named as gene names and 
  ##                 columns named as PC numbers. Potentially gained from Seurat Loadings() function
  
  ######### Varimax rotation
  initial_data <- t(gene_exp_matrix[rownames(gene_exp_matrix) %in% rownames(loading_matrix),])
  
  ## apply varimax rotation on the loadings
  varimax_res <- varimax(loading_matrix)
  rotatedLoadings <- varimax_res$loadings
  ## calculating the PC scores matrix
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  #scores          <- scale(initial_data) %*% invLoadings ## this second scaling is not necessary
  scores          <- initial_data %*% invLoadings ## this second scaling is not necessary
  ## compacting the rotated loading and score matrices in a list
  rotated_data <- list(rotLoadings=rotatedLoadings, rotScores = scores)
  return(rotated_data)
}




getDimRed_df <- function(seurat_obj, top_pc=1:30,
                         reduction='pca',
                         attributes=c('cell_type','sample','cluster')) {
  
  ## Extracting the top PCs for visualization and adding UMI attributes to the resulting dataframe
  ## INPUT:
  ###
  
  dimRed_df = data.frame(Embeddings(seurat_obj,reduction = reduction)[,top_pc])
  if(length(top_pc)==1) {
    colnames(dimRed_df)[1] = 'embedding_value'
  }
  
  if('cell_type' %in% attributes){
    dimRed_df$cell_type=seurat_obj$cell_type
  }
  if('sample' %in% attributes){
    dimRed_df$sample=seurat_obj$sample_name
  }
  if('cluster' %in% attributes){
    dimRed_df$cluster = Idents(seurat_obj)
  }
  return(dimRed_df)
}

