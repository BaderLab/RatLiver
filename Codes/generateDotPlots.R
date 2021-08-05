source('Codes/Functions.R')
Initialize()
#merged_samples <- readRDS('~/RatLiver/Objects/merged_samples_newSamples.rds')
merged_samples <- readRDS('~/RatLiver/Objects/merged_samples_newSamples.rds')
subset(pbmc, downsample = 100)
merged_samples <-FindVariableFeatures(merged_samples)


# Single cell heatmap of feature expression
DoHeatmap(subset(merged_samples, downsample = 1000), 
          features = VariableFeatures(merged_samples), size = 3) +theme(axis.title.y = element_blank())

#############  general annotation heatmaps #############

general_markers <- c('CYP3A7', 'CYP2A7', 'CYP2A6',
                    'SCD', 'HMGCS1', 'ACSS2', 'TM7SF2', 
                    'SEC16B', 'SLBP', 'RND3', 'PCK1', 
                    'BCHE', 'G6PC', 'GHR', 'ALDH6A1', 
                    'RPP25L', 'HSD11B1', 'HAMP', 'GHR',
                    'HPR', 'GSTA2', 'AKR1C1', 'MASP2' ,
                    'MGP', 'SPARCL1', 'TM4SF1', 'CLEC14A', 
                    'CCL14', 'CLEC1B', 'FCN2', 'S100A13', 
                    'RAMP3', 'INMT', 'DNASE1L3', 'LIFR',
                    'KRT7', 'KRT19', 'SOX9', 'EPCAM', 
                    'ACTA2', 'COL1A1', 'RBP1', 
                    'S100A8', 'LYZ', 'S100A9', 'HLA-DPB1', 
                    'CD5L', 'MARCO', 'VSIG4', 
                    'CD2', 'CD3D', 'TRAC', 'GZMK', 
                    'GNLY', 'PTGDS', 'GZMB', 'TRDC', 
                    'STMN1', 'HMGB2', 'TYMS', 
                    'MS4A1', 'LTB', 'CD37', 'CD79B', 
                    'IGLC2', 'IGHG1', 'IGKC', 
                    'CD7', 'KLRB1', 'NKG7', 
                    'HBB', 'CA1', 'ALAS2')


general_markers_df <- .getMapped_hs2model_df(ensembl, candidateGenes=general_markers, model_animal_name)
list_of_genes = c(general_markers_df$rnorvegicus_homolog_associated_gene_name)
is_unempty <- list_of_genes != '' & !is.na(list_of_genes)

data.frame(rat=list_of_genes[is_unempty],human=general_markers_df$symbol[is_unempty])
DotPlot(merged_samples, features = unique(list_of_genes[is_unempty])) + RotatedAxis()
DoHeatmap(merged_samples, features = unique(list_of_genes[is_unempty]), size = 3)



############# Checking the markers in the human liver paper - Supp info ############# 
Immune_marker = 'Cd68'
Inf_Mac = c('HLA-DRB5', 'HLA-DQA1', 'RP11-1143G9.4', 'CXCL8', 'VCAN', 'S100A12', 'MNDA', 'FCN1', 'IL18',
            'IL1R2', 'EVI2A', 'S100A9', 'S100A8', 'LYZ', 'C1QC', 'HLA-DRB1', 'HLA-DRA1', 'CD74', 'HLA-DRA', 'S100A6', 'TYROBP')

Non_Inf_Mac = c('CTSB', 'AIF1', 'CD163', 'VCAM1', 'CD5L', 'MARCO', 'MS4A7', 'LIPA', 'HMOX1', 'MAF', 
                'VSIG4', 'KLF4', 'C5AR1', 'CPVL', 'RAB31','VMO1', 'TMIGD3', 'LILRB5', 'TTYH3', 'CCDC88A', 'SLC31A2')


Antibody_secreting_B_cells = c('IGKC', 'SSR4', 'IGHG3', 'IGHA2', 'IGHG4', 'IGHG2', 'CD79A', 'RGS1', 'IL16', 'IGHGP', 'TMEM156', 
                               'FCRL5', 'PKHD1L1','IGLV3-1', 'IGLC7', 'JCHAIN', 'MZB1', 'IGHM', 'IGHG1', 'IGLC3', 'IGHA1', 'IGLC2')
Mature_B_cells = c('CD74', 'CD52', 'CD37', 'IRF8', 'TNFRSF13C', 'BANK1', 'VPREB3', 'IGHD', 'LINC00926', 'SELL',
                   'P2RX5', 'STAG3', 'RP5-887A10.1', 'RP11-693J15.5', 'TCL1A', 'CD22', 'CD79B', 'MS4A1', 'BIRC3', 'LTB', 'HLA-DQB1')
NK_like_cells = c('NKG7', 'XCL2', 'KLRC1', 'IL2RB', 'CD160', 'XCL1', 'CLIC3', 'TRDC', 'TXK', 'TMIGD2', 'CD69', 'GZMK', 'DUSP2',
                  'ALOX5AP', 'CMC1', 'KLRB1', 'CCL3', 'KLRF1', 'KLRD1', 'CD7')


gd_T_cells_1 <- c('NKG7', 'GZMA', 'CCL5', 'FGFBP2', 'GZMB', 'SPON2', 'HOPX', 'PTGDS', 'ADGRG1', 'S100B', 'CLIC3',
                'TRDC', 'GNLY', 'KLRF1', 'KLRD1', 'CD7', 'PRF1', 'CST7', 'IFITM1', 'CD247', 'CTSW')
gd_T_cells_2 <- c('GNLY', 'STMN1', 'TRDC', 'KIAA0101', 'TROAP', 'CENPA', 'ASPM', 'AURKB', 'CENPF', 'NUSAP1', 'MKI67',
                  'TOP2A', 'BIRC5', 'UBE2C', 'TYMS', 'H2AFX', 'PCNA', 'TUBB', 'HMGB2', 'TUBA1B')
ab_T_cells <- c('CD52', 'LTB', 'CD69', 'GZMK', 'TRBC1', 'TRBC2', 'DUSP2', 'ACAP1', 'CD3D', 'TRAC', 'CD2', 'CD3E', 
                'IL7R', 'PYHIN1', 'CCL4L2', 'EMB', 'TIGIT', 'IFNG', 'JUNB', 'IL32', 'CCL5', 'GZMA', 'HCST')

Hepatic_stellate_cells <- c('SPARC', 'BGN', 'TPM2', 'ACTA2', 'DCN', 'COL1A1', 'MYL9', 'COL3A1', 'COL1A2', 
                            'RBP1', 'CTGF', 'MEG3', 'HGF', 'CCL2', 'IGFBP3', 'IGFBP6', 'CYR61', 'OLFML3',
                            'COLEC11', 'IGFBP7', 'TAGLN')
Cholangiocytes <- c('DEFB1', 'CFTR', 'HNF1B', 'TFF3', 'LGALS2', 'SH3YL1', 'TACSTD2', 'KRT19', 'SCGB3A1', 
                    'CXCL1', 'CLDN4', 'CXCL6', 'TFF2', 'SOX9', 'EPCAM', 'MUC5B', 'TFF1', 'LGALS4', 'KRT7', 
                    'LCN2', 'ELF3', 'PIGR', 'CD24', 'SPP1', 'FXYD2')

Endothelial_cells <- c('TAGLN', 'FCGR2B', 'CLEC1B','FCN2', 'CCL14', 'DNASE1L3', 'PRSS23', 'LIFR', 'RAMP3', 
                       'S100A13', 'C7', 'CTGF', 'MGP', 'SPARCL1', 'CLEC14A', 'ADIRF', 'ID3', 'TM4SF1', 'PTGDS',
                       'PECAM1', 'SRPX', 'VWF', 'PCAT19', 'TSPAN7', 'F8', 'HYAL2', 'CALCRL', 'RAMP2', 'FLT1', 
                       'GJA4', 'NOSTRIN', 'STAB2', 'CD34', 'GATA4', 'CD36', 'LYVE1')

marker_ensembl_df <- .getMapped_hs2model_df(ensembl, candidateGenes=Endothelial_cells, model_animal_name)

list_of_genes = c(marker_ensembl_df$rnorvegicus_homolog_associated_gene_name)
is_unempty <- list_of_genes != '' & !is.na(list_of_genes)
data.frame(rat=list_of_genes[is_unempty],human=marker_ensembl_df$symbol[is_unempty])
DotPlot(merged_samples, features = unique(list_of_genes[is_unempty])) + RotatedAxis()


############# DotPlots for the top markers of each cluster #############
Cluster_markers <- readRDS('Results/new_samples/Cluster_markers.rds')
Cluster_markers_sorted <- lapply(Cluster_markers, function(x) x[order(x$p_val_adj, decreasing = F),])

i = 11
features <- rownames(Cluster_markers_sorted[[i]])[1:100]
Rp_genes_index <- grep(pattern = 'Rp', x = features)
Mt_genes_index <- grep(pattern = 'Mt-', x = features)

if(length(features[-c(Rp_genes_index, Mt_genes_index)])>0) features = features[-c(Rp_genes_index, Mt_genes_index)]
DotPlot(merged_samples, features = features) + 
  RotatedAxis() + ggtitle(names(Cluster_markers_sorted)[i])

features[grep(pattern = 'Rp', x = features)]


## subcluster each of the clusters mentioned in the slides
## based on the initial assumptions on the subcluster identity >> check the expression of a set of above markers using dotplots
## VAF project >> load the matrix and check the distribution

new_data_scCLustViz_object_Immune <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved_ImmuneSub.RData"
load(new_data_scCLustViz_object_Immune)
sCVdata_list$RNA_snn_res.1@Clusters
merged_samples <- your_scRNAseq_data_object

mac_markers <- read.csv('~/RatLiver/Results/new_samples/Mac_markers.csv')
table(mac_markers$Annotation)
mac_markers$Marker[!mac_markers$Marker %in% rownames(merged_samples)]

merged_samples$orig.ident <- as.character(sCVdata_list$RNA_snn_res.1@Clusters)
Idents(merged_samples) <- merged_samples$orig.ident
DotPlot(merged_samples, features = mac_markers$Marker[mac_markers$Marker %in% rownames(merged_samples)]) + 
  RotatedAxis() #+ ggtitle(names(Cluster_markers_sorted)[i])

