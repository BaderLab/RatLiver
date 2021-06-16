source('Codes/Functions.R')
Initialize()

### make expression input file (rows:genes, columns: samples(clusters))
merged_samples = readRDS('Objects/merged_samples_oldSamples_mt40_lib1500_MTremoved.rds')
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"

new_data_scCLustViz_object <- "Results/new_samples/scClustVizObj/for_scClustViz_newSamples_MTremoved.RData"

load(new_data_scCLustViz_object)
merged_samples <- your_scRNAseq_data_object
Idents(merged_samples) = as.character(sCVdata_list$res.0.6@Clusters)


Idents(merged_samples) = as.character(merged_samples$res.0.6)

cluster_names = as.character(Idents(merged_samples) )
cluster_names_types = unique(cluster_names)
mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)

### dividing the expression matrix based on the clusters
cluster_expression <- sapply(1:length(cluster_names_types), function(i){
  a_cluster_name = cluster_names_types[i]
  merged_samples[,cluster_names == a_cluster_name]
}, simplify = F)

names(cluster_expression) = cluster_names_types
lapply(cluster_expression, head)
names(cluster_expression) = paste0('cluster_', names(cluster_expression))
cluster_average_exp <- lapply(cluster_expression, function(x){
  ## calculate the average expression of each gene in each cluster
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})

lapply(cluster_average_exp, dim)
cluster_average_exp_df = do.call(cbind,cluster_average_exp)
colnames(cluster_average_exp_df) = names(cluster_average_exp)
head(cluster_average_exp_df)

row.names(cluster_average_exp_df) <- make.unique(df$V2)

dir_pathway <- 'Results/old_samples/pathway_analysis/'
dir_pathway <- 'Results/new_samples/pathway_analysis/'

dir.create(dir_pathway)
write.csv(cluster_average_exp_df, 
          paste0(dir_pathway, 'cluster_average_exp_merged_new_samples.txt')) # cluster_average_exp_merged_old_samples.txt





## make the annotation input file (GeneType file indicates which gene is protein coding)
dir = '~/Desktop/gsva_inputs/'
dir.create(dir)
dir = 'Codes/pathway_analysis/gsva/'

gft <- rtracklayer::import("~/Desktop/Rattus_norvegicus.Rnor_6.0.100.gtf")
gtf_df = as.data.frame(gft)
head(gtf_df)

#gtf_df_unique <- gtf_df[!duplicated(gtf_df$gene_id),]
dim(gtf_df_unique)
gtf_df_unique <- gtf_df
gtf_df_unique <- gtf_df_unique[,c('gene_id', 'gene_name', 'gene_biotype')]
head(gtf_df_unique)

write.table(gtf_df_unique, paste0(dir,'GeneType.txt')) 

