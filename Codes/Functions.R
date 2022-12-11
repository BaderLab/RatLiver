## This script includes functions needed in the analysis 
## and need to be loaded in the beginning of scripts

# function to check to see if packages are installed. 
#  Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


## loading required packages    
Initialize <- function(){
  options(stringsAsFactors = FALSE)
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  options(stringsAsFactors = FALSE)
  listOfPackages <- c('BiocManager', 'Seurat', 'viridis', 'org.Rn.eg.db', 'stringr',  'biomaRt', 
                      'ggplot2', 'scClustViz', 'org.Hs.eg.db', 'AnnotationDbi', 'presto' ,
                      'scran', 'Matrix', 'devtools', 'AUCell', 'GSEABase','GSVA', 'fgsea','limma', 'SCINA',
                      'DelayedArray', 'DelayedMatrixStats','monocle', 'multtest', 'gprofiler2',
                      'stringr', 'harmony', 'SoupX', 'pheatmap') #, 'celda' >> installation issue
  ipak(listOfPackages)
  
}

library('Seurat')
library('org.Rn.eg.db')
library('scClustViz')


colorPalatte <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "brown"
)

getHead <- function(dataframe){print(dataframe[1:5, 1:5])}

getUnemptyList <- function(chrList){ chrList[!is.na(chrList) & chrList != '' ]}
getUnemptyList_bool <- function(chrList){ !is.na(chrList) & chrList != '' }


get_cluster_object <- function(seur, 
                               max_seurat_resolution,
                               FDRthresh,
                               min_num_DE,
                               seurat_resolution,
                               seurat_resolution_jump,
                               DE_bw_clust){
  sCVdata_list <- list()
  
  while(DE_bw_clust) {
    if (seurat_resolution >= max_seurat_resolution) { break }
    seurat_resolution <- seurat_resolution + seurat_resolution_jump 
    # ^ iteratively incrementing resolution parameter 
    
    seur <- FindClusters(seur,resolution=seurat_resolution,verbose=F)
    
    message(" ")
    message("------------------------------------------------------")
    message(paste0("--------  res.",seurat_resolution," with ",
                   length(levels(Idents(seur)))," clusters --------"))
    message("------------------------------------------------------")
    
    if (length(levels(Idents(seur))) <= 1) { 
      message("Only one cluster found, skipping analysis.")
      next 
    } 
    # ^ Only one cluster was found, need to bump up the resolution!
    
    if (length(sCVdata_list) >= 1) {
      temp_cl <- length(levels(Clusters(sCVdata_list[[length(sCVdata_list)]])))
      if (temp_cl == length(levels(Idents(seur)))) { 
        temp_cli <- length(levels(interaction(
          Clusters(sCVdata_list[[length(sCVdata_list)]]),
          Idents(seur),
          drop=T
        )))
        if (temp_cli == length(levels(Idents(seur)))) { 
          message("Clusters unchanged from previous, skipping analysis.")
          next 
        }
      }
    }
    
    curr_sCVdata <- CalcSCV(
      inD=seur,
      assayType="RNA",
      assaySlot="counts",
      cl=Idents(seur), 
      # ^ your most recent clustering results get stored in the Seurat "ident" slot
      exponent=NA, 
      # ^ going to use the corrected counts from SCTransform
      pseudocount=NA,
      DRthresh=0.1,
      DRforClust="harmony", #pca
      calcSil=T,
      calcDEvsRest=T,
      calcDEcombn=T
    )
    
    DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
    # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
    message(paste("Number of DE genes between nearest neighbours:",min(DE_bw_NN)))
    
    if (min(DE_bw_NN) < min_num_DE) { DE_bw_clust <- FALSE }
    # ^ If no DE genes between nearest neighbours, don't loop again.
    
    sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
  }
  return(sCVdata_list)
}




Plot.tsne.gene.expr <- function(tsne.gene.df, title='',subtitle=''){
  colnames(tsne.gene.df)[3] <- 'expression'
  ggplot(tsne.gene.df, aes(x=tSNE_1, y=tSNE_2, color=expression))+
    geom_point(alpha=0.7)+theme_classic()+ggtitle(title,subtitle = subtitle) + 
    theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, color = "blue"))+ 
    xlab('tSNE1')+ylab('tSNE2')+
    scale_color_viridis(direction = -1,option = "plasma") #limits=c(0,3),oob=scales::squish
}


getTsneDF <- function(seur){
  listOfClusters <- as.character(Idents(seur))
  DataFrame <- as.data.frame(cbind(getEmb(seur, 'tsne'),clusters=listOfClusters))
  DataFrame$tSNE_1 <- as.numeric(DataFrame$tSNE_1)
  DataFrame$tSNE_2 <- as.numeric(DataFrame$tSNE_2)
  return(DataFrame)
}

getUmapDF <- function(seur){
  listOfClusters <- as.character(Idents(seur))
  DataFrame <- as.data.frame(cbind(getEmb(seur, 'umap'),clusters=listOfClusters))
  DataFrame$UMAP_1 <- as.numeric(DataFrame$UMAP_1)
  DataFrame$UMAP_2 <- as.numeric(DataFrame$UMAP_2)
  return(DataFrame)
}

getPcaDF <- function(seur){
  listOfClusters <- as.character(Idents(seur))
  DataFrame <- as.data.frame(cbind(getEmb(seur, 'pca')[,1:2],clusters=listOfClusters))
  DataFrame$PC_1 <- as.numeric(DataFrame$PC_1)
  DataFrame$PC_2 <- as.numeric(DataFrame$PC_2)
  return(DataFrame) 
}

detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


get_manual_labels <- function(){
  
  ### setting the cluster names manually
  clusters_rat_DA_01 <- data.frame(cluster=paste0('cluster_', c(0:11)))
  clusters_rat_DA_01$cell_type <- c(rep('Hepatocyte',4),'KC', 'LSEC','Stellate_cell','T_cell',
                                    'Immune_cell', 'T_cell','LSEC','DC_cell')
  
  clusters_rat_DA_02 <- data.frame(cluster=paste0('cluster_', c(0:18)))
  clusters_rat_DA_02$cell_type <- c(rep('Hepatocyte',6),'LSEC','Hepatocyte','KC', rep('Hepatocyte',2),
                                    'Stellate_cell','Immune_cell','T_cell','T_cell','Immune_cell',
                                    'DC_cell','Prolif_cell','Erythrocyte')
  
  clusters_rat_Lew_01 <- data.frame(cluster=paste0('cluster_', c(0:15)))
  clusters_rat_Lew_01$cell_type <- c(rep('Hepatocyte',4),'LSEC','Hepatocyte','KC','Stellate_cell',
                                     'KC','T_cell', 'Hepatocyte', 'Immune_cell','T_cell', 'Hepatocyte',
                                     'Immune_cell','DC_cell')
  
  clusters_rat_Lew_02 <- data.frame(cluster=paste0('cluster_', c(0:14)))
  clusters_rat_Lew_02$cell_type <- c(rep('Hepatocyte',4), 'KC', 'Hepatocyte','LSEC','Hepatocyte', 
                                     'Stellate_cell','KC','T_cell', 'DC_cell', 'T_cell', 'Immune_cell', 
                                     'Stellate_cell')
  
  
  clusters_rat_DA_M_10WK_003 <- data.frame(cluster=paste0('cluster_', c(0:9)))
  clusters_rat_DA_M_10WK_003$cell_type = c(rep('Hepatocyte',2), 'LSEC','Hepatocyte' ,'KC','Stellate_cell',
                                           'Hepatocyte','Immune_cell','T_cell','DC_cell' )
  
  
  cluster_cell_type_df <- list(clusters_rat_DA_M_10WK_003, clusters_rat_DA_01, 
                               clusters_rat_Lew_01, clusters_rat_Lew_02, clusters_rat_DA_02)
  names(cluster_cell_type_df) = c('rat_DA_M_10WK_003','rat_DA_01', 'rat_Lew_01', 'rat_Lew_02', 'rat_DA_02')
  
  return(cluster_cell_type_df)
}


colorPalatte <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "brown"
)


get_vector_of_lables <- function(Cell_type_assigned, seur){
  ifelse(colnames(seur) %in% Cell_type_assigned$STELLATE, 'STELLATE',
         ifelse(colnames(seur) %in% Cell_type_assigned$B_CELLS, 'B_CELLS', 
                ifelse(colnames(seur) %in% Cell_type_assigned$M1_KC, 'M1_KC', 
                       ifelse(colnames(seur) %in% Cell_type_assigned$M2_KC, 'M2_KC', 
                              ifelse(colnames(seur) %in% Cell_type_assigned$CYTOTOXIC_CELLS, 'CYTOTOXIC_CELLS', 
                                     ifelse(colnames(seur) %in% Cell_type_assigned$HEMATOPOIETIC, 'HEMATOPOIETIC', 
                                            ifelse(colnames(seur) %in% Cell_type_assigned$T_CELLS, 'T_CELLS',
                                                   ifelse(colnames(seur) %in% Cell_type_assigned$LSECS, 'LSECS', 
                                                          ifelse(colnames(seur) %in% Cell_type_assigned$CHOLANGIOCYTES, 'CHOLANGIOCYTES', 
                                                                 ifelse(colnames(seur) %in% Cell_type_assigned$HEPATOCYTE, 'HEPATOCYTE', 'unknown'
                                                                 ))))))))))
  
}


#breaks=c(0,0.5,1),labels=c("Minimum",0.5,"Maximum"),
# limits=c(0,1)


#plot_tsne(cell_coord=getEmb(seur,"tsne"),
#          md=gene.expr,
#          md_title=GENE_NAME,
#          md_log=F)
