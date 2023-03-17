library(Seurat)

path_to_scClustViz_obj <- "~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC_scCLustViz_object_res.0.6.RData"
load(path_to_scClustViz_obj)

resolutions = 0.6
merged_samples <- FindClusters(merged_samples, resolution = resolutions, verbose = FALSE)
Idents(merged_samples) = merged_samples$SCT_snn_res.0.6


num_top_genes = 20
path_to_markers_list = '~/rat_sham_sn_data/standardQC_results/cluster_markers_res0.6/cluster_8.txt'
markers_df <- read.csv(path_to_markers_list)
head(markers_df)
gene_names = unique(markers_df$X[1:num_top_genes])

DotPlot(merged_samples, features = gene_names) + RotatedAxis()
