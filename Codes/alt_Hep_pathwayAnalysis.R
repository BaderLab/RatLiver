#GSVA for samples, please keep the .gmt file in the same folder
require(GSVA)
require(parallel)
require(GSA)
require(qvalue)
library('qusage')

#Input data with gene symbol as the first column (name ID) and sample names as the first row
# f1 <- readline(prompt="Please enter the full name (including .txt, .csv, etc)
#                of the gene expression file (gene ID 1st column, sample name 1st row): ")
# datam <- read.table( f1, header = TRUE, sep = "\t", quote="\"",  stringsAsFactors = FALSE)

f1 = 'Results/rat_old_cluster_average_exp_all.rds'
datam <- readRDS(f1)
datam <- cbind(genes=rownames(datam),datam[,-ncol(datam)])
head(datam)

#f3 <- readline(prompt="Is the Data: 1. normalized microarray/RNA-seq;
#             or 2. RNA-seq raw counts (Input 1 or 2): ")
f3='1'
if (f3=="2") {rseq="Poisson"} else {rseq="Gaussian"}

f4 <- 'Results/resources/GeneType.txt'
genetype <- read.table(f4, header = T)
proteingene <- genetype[which(genetype[,3]=="protein_coding"),1:2]
#remove all duplicated genes
proteingene <- proteingene[!duplicated(proteingene[,2]),]

#Only protein coding and unique genes used
datau <- datam[which(datam[,1]%in%proteingene[,2]),]
datau <- datau[!duplicated(datau[,1]),]
rownames(datau) <- proteingene[match(datau[,1],proteingene[,2]),2]
datau <- as.matrix(datau[,2:ncol(datau)])
head(datau)
dim(datau)

cell_type_names = unlist(lapply(str_split(colnames(datau), pattern = ' '),'[[', 1))
hep_clusters <- which(cell_type_names=='Hep')

### selecting the Hep clusters
datau <- datau[,hep_clusters] #set2_heps
head(datau)
Hep_0_15_1 = rowMeans(datau[,paste0('Hep (', c(0, 15, 1), ')')])
Hep_6_12_4 = rowMeans(datau[,paste0('Hep (', c(6, 12, 4), ')')])
Hep_16_2_8 = rowMeans(datau[,paste0('Hep (', c(16, 2, 8), ')')])

datau = cbind(Hep_0_15_1, Hep_6_12_4, Hep_16_2_8)
head(datau)

#Input the file with subtype information for each sample, first column contains sample name, second column records group name
f2 <- list.files('Results/resources/', pattern = ".gmt", full.names = T)

gmt1 <- GSA.read.gmt(f2)
gmt2 <- list()
for (i in 1:length(gmt1$genesets)){
  gene_list_i <- unlist(gmt1$genesets[i])
  description_gene_set_i <- gmt1$geneset.descriptions[i]
  gmt2[[description_gene_set_i]] <- gene_list_i[gene_list_i!=""]
}

table(data.frame(str_split_fixed(gmt1$geneset.names,':',2))[,1])

#GSVA and output
result <- gsva(datau, gmt2, min.sz=10, max.sz=200, mx.diff=TRUE,
               verbose=T, kcdf=rseq, parallel.sz=0)# ,,method='ssgsea'
result1 <- as.matrix(result)
head(result1)

result2 <- NULL
for (i in 1:nrow(result1)){
  # go through rows of results matrix and find the correspoing GO-id for the GO-description
  # gene_set_id <- unlist(strsplit(gmt1$geneset.names[match(rownames(result1)[i],
  #                                               gmt1$geneset.descriptions)],"%"))[3]
  gene_set_id <- gmt1$geneset.names[gmt1$geneset.descriptions == rownames(result1)[i]]
  gene_names_pasted <- paste(unlist(gmt2[[match(rownames(result1)[i],names(gmt2))]]),collapse = ",")
  result2 <- rbind(result2,c(gene_set_id,gene_names_pasted))
}

dir = "Results/old_samples/pathway_analysis/alt_Hep_pathway/"
#dir = "~/gsva_inputs/GSVA_Gprofile/"

dir.create(dir)
for (j in 1:ncol(result1)){
  
  p_value_list <- pnorm(-abs(scale(result1[,j])))
  q_value_list <- qvalue(p_value_list,lambda=seq(0.05,0.45,0.01))$lfdr
  gene_set_ids <- result2[,1]
  gene_set_descriptions <- rownames(result1)
  gene_names_pasted_list <- result2[,2]
  
  final_df <- cbind(gene_set_ids, gene_set_descriptions,
                    p_value_list, q_value_list, sign(result1[,j]), gene_names_pasted_list)
  colnames(final_df)<-c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")
  
  write.table(final_df,
              paste0(dir,colnames(result1)[j],"_gProfile.txt"),
              sep = "\t",row.names = F,quote = F)
  
  write.table(final_df[final_df[,5]>0,],
              paste0(dir,colnames(result1)[j],"_gProfilePos.txt"),
              sep = "\t",row.names = F,quote = F)
}

#write.csv(result1,file=paste0(f1,"_ssgsea_protein_result.csv"))
write.csv(result1,file=paste0(f1,'_GSVA_protein_result.csv'))


####### performing differential analysis between clusters 
# periportal: (8,2,16)  pericentral: (4, 12, 6, 1, 15, 0)


old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
merged_samples <- your_scRNAseq_data_object
merged_samples$clusters <- as.character(sCVdata_list$res.0.6@Clusters)
merged_samples <- merged_samples[,merged_samples$clusters %in% c(4, 12, 6, 1, 15, 0, 8,2,16)]
dim(merged_samples)
dim(your_scRNAseq_data_object)

table(merged_samples$clusters)
Idents(merged_samples) <- ifelse(merged_samples$clusters %in% c(8,2,16), 'pp', 'cv')
table(Idents(merged_samples))

categories = levels(merged_samples)
pp_cv_markers <- sapply(1:length(categories), 
                                  function(i) FindMarkers(merged_samples, 
                                                          ident.1=categories[i],
                                                          logfc.threshold = 0,
                                                          min.pct=0,
                                                          min.cells.feature = 1,
                                                          min.cells.group = 1), 
                                  simplify = FALSE)

saveRDS(pp_cv_markers, 'Results/pp_cv_markers.rds')
names(pp_cv_markers) = categories

pp_cv_markers_ord = lapply(pp_cv_markers, function(x){
  x$ranking_score=-log10(x$p_val_adj+.Machine$double.xmin)*sign(x$avg_log2FC)
  x = x[order(x$ranking_score, decreasing = T),]
  return(x)
})

pp_cv_markers_ord = lapply(pp_cv_markers, function(x){
  x[order(x$avg_log2FC, decreasing = T),]
})
lapply(pp_cv_markers_ord, head, 30)

dir.create('~/RatLiver/Results/pp_cv_set1Hep_pathway/')
for(i in 1:2){
  df = data.frame(genes=rownames(pp_cv_markers_ord[[i]]), score=pp_cv_markers_ord[[i]]$avg_log2FC)
  
  write.table(df, 
              paste0('~/RatLiver/Results/pp_cv_set1Hep_pathway/',categories[i],'.rnk'), 
              #paste0('~/RatLiver/Results/old_samples/rotated_loadings/rot_PC',pc_index,'_loadings.rnk'), 
              col.names = F, row.names = F, quote = F, sep = '\t')
  
  
}


source('Codes/Functions.R')
Initialize()
gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'



####### Running GSEA on the markers list ########
rnk_file_path = '~/RatLiver/Results/pp_cv_set1Hep_pathway/'
rnk_files <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = T)
working_dir = paste0(rnk_file_path, '/gsea_results/')
dir.create(working_dir)

for(i in 1:length(rnk_files)){
  rnk_file = rnk_files[i]
  print(rnk_file)
  analysis_name = gsub('.rnk', '',list.files(rnk_file_path,pattern = '*.rnk')[i])
  working_dir = working_dir
  
  GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                       rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                       analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , #set min=15
                       working_dir, " > gsea_output.txt")
  system(GSEA_command)
}

