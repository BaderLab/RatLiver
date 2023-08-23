# script JL-GSVAGprofileproteingeneUniquePath.R
# needs you put few files in a folder. 1: the Cluster Average matrix by genes;
# 2: your gmt file;
# 3: Gene type annotation.
# I attached the human version of Gene annotation here so you can generate one for rat.
# The script will filter out non-protein genes first to speed up processing.
# I also notice that different gsva version may lead to errors.
# If your version of GSVA is generating errors, just delete the "kcdf" part of the code:
#
# result<-gsva(datau, gmt2, min.sz=10, max.sz=200, mx.diff=TRUE, verbose=F, kcdf=rseq, parallel.sz=0)
#
# There are also old version of gsva will give error for:
#   result1<-as.matrix(result) because result will be in different format.
#
# I format the output of the analysis into G:profiler output so we can use Cytoscape to visualize the results.
# As for p-values, you need another script JL-TCGA-GSVA.R where bootstrap will be performed.

#GSVA for samples, please keep the .gmt file in the same folder

# I have calculated P-values and FDR without bootstrap, but it is by distribution of the enrichment scores.
# The assumption is that most of the pathways will not be enriched and the ones with high NES will be
# considered significant. This is usually pretty good but bootstrap is a better approach, it randomized
# the gene expression data and perform multiple GSVA (>100) to make sure that for that pathway under
# the same data the NES calculated is significant. It is up to you to perform the bootstrap or not.

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


f1 = '~/rat_sham_sn_data/standardQC_results/sham_sn_merged_annot_standardQC_averageCluterExp.csv'
datam <- read.csv(f1, header = T)
head(datam)
#f3 <- readline(prompt="Is the Data: 1. normalized microarray/RNA-seq;
#             or 2. RNA-seq raw counts (Input 1 or 2): ")
f3='1'
if (f3=="2") {rseq="Poisson"} else {rseq="Gaussian"}


f4 <- '~/RatLiver/Results/resources/GeneType.txt'
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

#Input the file with subtype information for each sample, first column contains sample name, second column records group name
f2 <- list.files('~/RatLiver/Results/resources/', pattern = ".gmt", full.names = T)

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

dir = '~/rat_sham_sn_data/standardQC_results/gsva_result/'
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


