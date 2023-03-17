source('Codes/Functions.R')
Initialize()
################################
# Running  pathway analysis in the nuc seq data hepatocyte populations
################################

merged_samples = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC.rds')
resolutions = 0.6
merged_samples <- FindClusters(merged_samples, resolution = resolutions, verbose = FALSE)

table(merged_samples$annot_IM_g,merged_samples$SCT_snn_res.0.6)
table(merged_samples$SCT_snn_res.0.6[merged_samples$annot_IM_g == 'Hep'])
nucseq_heps = paste0('cluster_', as.character(c(0:4, 6, 9:11)))


################## importing the list of markers for each population
Cluster_markers_final <- readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_cluster_markers_res0.6.rds')
length(Cluster_markers_final)
names(Cluster_markers_final)

#### subsetting the cluster markers to include the Hepatocytes
Cluster_markers_final <- Cluster_markers_final[names(Cluster_markers_final) %in% nucseq_heps]
names(Cluster_markers_final)

#### saving markers #####
dir.create(paste0('~/rat_sham_sn_data/standardQC_results/cluster_markers_res0.6_rnk_files/'))

for(i in 1:length(Cluster_markers_final)){
  
  df <- data.frame(Cluster_markers_final[[i]])
  
  df$score = -log10(df$p_val_adj+1e-50)*df$avg_log2FC
  df_ord = df[order(df$score, decreasing = T),]
  #df_ord = df[order(df$avg_log2FC, decreasing = T),]
  
  df_ord = data.frame(gene=rownames(df_ord), score=df_ord$score)
  head(df_ord,20)
  
  file_name <- names(Cluster_markers_final)[i]
  
  write.table(df_ord, 
              paste0('~/rat_sham_sn_data/standardQC_results/cluster_markers_res0.6_rnk_files/', file_name,'.rnk'),
              quote = F, sep = '\t',row.names = FALSE, col.names = FALSE)
}




##### Run GSEA on the DE list



####### Running GSEA on the markers list ########
gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'


rnk_file_path = paste0('~/rat_sham_sn_data/standardQC_results/cluster_markers_res0.6_rnk_files/')
rnk_files <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = T)
working_dir = paste0('~/rat_sham_sn_data/standardQC_results/cluster_markers_res0.6_rnk_files/gsea_results/')
dir.create(working_dir)


for(i in 2:length(rnk_files)){
  rnk_file = rnk_files[i]
  print(rnk_file)
  
  analysis_name = gsub('.rnk', '',list.files(rnk_file_path, pattern = '*.rnk')[i])
  working_dir = working_dir
  
  GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                       rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                       analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , #set min=15
                       working_dir, " > gsea_output.txt")
  system(GSEA_command)
  }



