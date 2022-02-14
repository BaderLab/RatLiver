### look into the expression pattern of CD68, CD11b/c, CD206, CD163
#### check the markers of that populaton vs all the other cells and vs other MAC populations

source('Codes/Functions.R')
############
old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)

#### the new immune subcluster data ####
merged_samples <- your_scRNAseq_data_object
merged_samples$cluster = as.character(sCVdata_list$res.0.6@Clusters)

merged_samples$sample_name = ifelse(merged_samples$orig.ident=='rat_DA_01_reseq', 'DA-1', 
                                    ifelse(merged_samples$orig.ident=='rat_DA_M_10WK_003', 'DA-2',
                                           ifelse(merged_samples$orig.ident=='rat_Lew_01', 'Lew-1', 'Lew-2')))
merged_samples$strain = sapply(str_split(colnames(merged_samples), '_'), '[[', 2)
#

#### check if removal of mt-genes would help with the annotation
df_umap <- data.frame(UMAP_1=getEmb(merged_samples, 'umap')[,1], 
                      UMAP_2=getEmb(merged_samples, 'umap')[,2], 
                      clusters=merged_samples$cluster,
                      sample=merged_samples$sample_name,
                      strain=merged_samples$strain,
                      Ptprc=GetAssayData(merged_samples)['Ptprc',], 
                      umi=colnames(merged_samples))

df_umap$Cd14 <-  GetAssayData(merged_samples)['Cd14',]
##### adding the final annotations to the dataframe
set1_info <- read.csv('figure_panel/set-1-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_info$clusters = as.character(set1_info$clusters)
df_umap = merge(df_umap, set1_info[,1:3], by.x='clusters', by.y='clusters', all.x=T,order=F)
#### re-ordering the rows based on UMIs of the initial data
df_umap <- df_umap[match(colnames(merged_samples),df_umap$umi),]
df_umap$umi == colnames(merged_samples)


# Mrc1 == Cd206
# Itgam = Cd11b --> not in the dataset
# Itgax = Cd11c
sai_genes = c('Cd68', 'Mrc1', 'Itgax','Cd163', 'Itgal')
sai_genes %in% rownames(merged_samples)
rownames(merged_samples)[grep('Itg', rownames(merged_samples))]
sai_genes.df = t(GetAssayData(merged_samples)[rownames(merged_samples)%in%sai_genes,])

############## Varimax-related plots ##############
rot_data <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('Varimax_', 1:ncol(scores))
embedd_df_rotated <- data.frame(scores)

nb.cols <- length(unique(df_umap$label))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

pc_num = 15
rot_df <- data.frame(Varimax_1=embedd_df_rotated$Varimax_1,
                     emb_val=embedd_df_rotated[,pc_num],
                     cluster=as.character(merged_samples$cluster),
                     Library_size=merged_samples$nCount_RNA,
                     num_expressed_genes=merged_samples$nFeature_RNA,
                     strain=merged_samples$strain,
                     sample_name = merged_samples$sample_name,
                     label=factor(df_umap$label),
                     sai_genes.df) #, levels = as.character(set_1_cl_ord)

head(rot_df)
ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=label))+geom_point(alpha=0.7,size=2)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(name='cell-type',values = mycolors)+theme_classic()+
  theme(text = element_text(size=22))+#, legend.title = element_blank()
  ggtitle(paste0('Set-1 cells over Varimax 1 and ', pc_num))



tnfr_genes = rownames(merged_samples)[grep(rownames(merged_samples), pattern = 'Tnfr')]

i = 4
gene_name = tnfr_genes[i]
rot_df$a_gene <-  GetAssayData(merged_samples)[gene_name,]
ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=a_gene))+geom_point(alpha=0.7,size=1.3)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_viridis(direction = -1)+theme_classic()+
  theme(text = element_text(size=18),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=13,angle=90,color='black'),
        legend.title = element_blank(),legend.position="top")+
  guides(colour=guide_colourbar(barwidth=15))+
  ggtitle(paste0(gene_name))

rownames(merged_samples)[grep('Cd1', rownames(merged_samples))]

#### expression/score over UMAP
a_gene = 'Itgal' 
#df_umap$a_gene = abs(rot_df$emb_val)
df_umap$a_gene =  GetAssayData(merged_samples)[a_gene,]
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(alpha=0.5,size=1)+
  scale_color_viridis(direction = -1, option = 'inferno',name=a_gene)+theme_classic()+ #paste0("Varimax-", pc_num)
  theme(text = element_text(size=22), legend.title=element_text(size=17))


######### checking the relationship between Itgal and 
df_umap$Itgal =  GetAssayData(merged_samples)['Itgal',]
df_umap$Mrc1 =  GetAssayData(merged_samples)['Mrc1',]
df_umap2 = df_umap[df_umap$label %in% c('Non-Inflammatory Mac (10)', 'Non-Inflammatory Mac (5)'),]
ggplot(df_umap2, aes(x=Itgal, y=Mrc1))+stat_binhex()+scale_fill_gradient(low = "blue", high = "red")+
  geom_smooth(method = "lm", se = FALSE)
 

# scale_fill_gradient(low = "lightblue", high = "red")+stat_bin2d(bins = 50) 

####### Check the enrichment of TNFa-related genes in each strain
geneSets <- readRDS('~/geneSets_converted2rat.rds') ### CHEA genesets converted to rat genes
a_gene_set_name <- 'TNF' # 'HNF4A' # GATA1
query_genesets = names(geneSets)[grepl(pattern = a_gene_set_name, names(geneSets),ignore.case = T)]


gmtFile = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'
geneSets <- getGmt(gmtFile)
names(geneSets)
query_genesets = names(geneSets)[grepl(pattern = 'TNF[-]?A', names(geneSets),ignore.case = T)][-3]
my.matrix.sub = GetAssayData(merged_samples)[,df_umap$label %in% c('Non-Inflammatory Mac (10)','Non-Inflammatory Mac (5)')]

TNFa_genesets = geneIds(geneSets[names(geneSets) %in% query_genesets])




#### evaluating the genesets:
# https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html
# > Genes regulated by NF-kB in response to TNF [GeneID=7124].
MsigDB = read.table('MsigDB_halmark_TNF_sig_via_NFLB.txt', header = T)
candidateGenes= MsigDB$HALLMARK_TNFA_SIGNALING_VIA_NFKB
source('Codes/convert_human_to_ortholog_functions.R')
MsigDB.mus = .getMapped_hs2model_df(ensembl, candidateGenes, model_animal_name)
MsigDB.mus.genes = MsigDB.mus$rnorvegicus_homolog_associated_gene_name[MsigDB.mus$rnorvegicus_homolog_orthology_type=='ortholog_one2one']

TNFa_genesets$MsigDB.TNFa.manual = MsigDB.mus.genes



########## evaluating it for Kegg's TNF geneset
library("KEGGREST")
kegg_TNF = keggGet("rno04668")[[1]] # "TNF signaling pathway - Rattus norvegicus (rat)"
kegg_TNF$REFERENCE
kegg_TNF$GENE

#Get the list of numbers, gene symbols and gene description
names <- kegg_TNF$GENE
#Delete the gene number by deleting every other line
namesodd <-  names[seq(0,length(names),2)]
#Create a substring deleting everything after the ; on each line (this deletes the gene description).
namestrue <- gsub("\\;.*","",namesodd)
#export the vector as a csv

TNFa_genesets$Kegg.TNFa.manual = namestrue

############################ 
lps_tnf_df = read.csv('lps_tnf_genes.csv')
tnf_secretion.geneset = getUnemptyList(lps_tnf_df$rnorvegicus_homolog_associated_gene_name[lps_tnf_df$TNF=='+'])

TNFa_genesets$tnf_secretion = tnf_secretion.geneset


############################ CMap purturb tnf genes
petrurb_tns.rat = read.csv('tnf_CMap_purturb_genes.csv')
TNFa_genesets$tnf_purturb = getUnemptyList(petrurb_tns.rat$rnorvegicus_homolog_associated_gene_name)


############################
scores_2 <- data.frame(ScoreSignatures_UCell(my.matrix.sub, features=TNFa_genesets))
scores_2$strain = merged_samples$strain[df_umap$label %in% c('Non-Inflammatory Mac (10)','Non-Inflammatory Mac (5)')]
head(scores_2)

NUM_GENE_SETS = 9
colnames(scores_2)[1:NUM_GENE_SETS] = c('TNF_alpha (IOB)', 'Halmark_TNFA_signaling_via_NFKB (MSIGDB_C2)',
                            'TNF_alpha (NETPATH)','TNF_a_NFKB_signaling (WIKIPATHWAYS)',
                            'TNFa_mucus_production (WIKIPATHWAYS)', 'TNFa (MSIGDB_manual)', 'TNFa (KEGG_manual)',
                            'TNF secretion', 'TNF purturbed (CMap)')


names(TNFa_genesets) = colnames(scores_2)[1:NUM_GENE_SETS]
lapply(TNFa_genesets, length)
#### genes in each TNFa geneset which are not present in the dataset 
lapply(TNFa_genesets, function(x) x[!x%in%rownames(merged_samples)])
#### genes in each TNFa geneset which are included in the rat liver dataset 
TNFa_genesets = lapply(TNFa_genesets, function(x) x[x%in%rownames(merged_samples)])

######## calculating the scores for the whole dataset
my.matrix = GetAssayData(merged_samples)
scores_3 <- data.frame(ScoreSignatures_UCell(my.matrix, features=TNFa_genesets))
scores_3$strain = merged_samples$strain
colnames(scores_3)[1:NUM_GENE_SETS] = names(TNFa_genesets)

head(scores_3)


##### evaluating the strain differences between the strains

pdf('TNFa_genesetUcell_strainSpecific_3.pdf', width = 15,height = 10)
for( i in 1:NUM_GENE_SETS){
  print(i)
  data =  data.frame(PathwayScore=scores_2[[i]], 
                     strain=scores_2$strain)
  
  p1=ggplot(data, aes(y=PathwayScore, x=strain))+geom_boxplot(aes(fill=strain))+
    ylab(colnames(scores_2)[i])+theme_classic()+
    scale_fill_manual(values=c("#E69F00", "#009E73"))+ggtitle(colnames(scores_2)[i])+
    stat_summary(fun.data=data_summary)+stat_compare_means()#method = "t.test"
  
  

  
  geneset.df = data.frame(t(my.matrix[TNFa_genesets[[i]],]))
  geneset.df$strain = merged_samples$strain
  head(geneset.df)
  geneset.df.mean = data.frame(Lew=colSums(geneset.df[geneset.df$strain=='Lew',-ncol(geneset.df)])/sum(geneset.df$strain=='Lew'),
             DA=colSums(geneset.df[geneset.df$strain=='DA',-ncol(geneset.df)])/sum(geneset.df$strain=='DA'))
  
  head(geneset.df.mean)
  geneset.df.mean$dif = geneset.df.mean$Lew - geneset.df.mean$DA
  geneset.df.mean = geneset.df.mean[order(geneset.df.mean$dif, decreasing = T),]
  
  LEW_TNFa_genes = head(geneset.df.mean[order(geneset.df.mean$dif, decreasing = T),], 20)
  DA_TNFa_genes = head(geneset.df.mean[order(geneset.df.mean$dif, decreasing = F),], 20)
  
  head(LEW_TNFa_genes)
  head(DA_TNFa_genes)
  
  lim = 20
  if(nrow(LEW_TNFa_genes)<lim) lim = nrow(LEW_TNFa_genes)
  genes2Vis = c(rownames(LEW_TNFa_genes)[1:lim], rownames(DA_TNFa_genes)[1:lim])
  Idents(merged_samples) <- factor(merged_samples$strain)
  p4=DotPlot(merged_samples, features = unique(genes2Vis)) + RotatedAxis()+ xlab('Markers')+ylab('')+
    theme(axis.text.x = element_text(size=9)) + 
    ggtitle(paste0('strain-specific genes within the ', colnames(scores_2)[i]))

  
  
  
  #### vairmax plot
  rot_df$a_geneset <-  scores_3[[i]]
  p2=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=a_geneset))+geom_point(alpha=0.7,size=1.3)+
    theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
    scale_color_viridis(direction = -1, option='inferno')+theme_classic()+
    theme(text = element_text(size=9),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=9,angle=90,color='black'),
          legend.title = element_blank(),legend.position="top")+
    guides(colour=guide_colourbar(barwidth=15))+
    ggtitle(colnames(scores_2)[i])

  
  ##### UMAP plot
  df_umap$a_geneset <-  scores_3[[i]]
  p3=ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_geneset))+geom_point(alpha=0.5,size=1)+
    scale_color_viridis(direction = -1, option = 'inferno',name='Ucell TNF')+theme_classic()+ #paste0("Varimax-", pc_num)
    theme(text = element_text(size=9), legend.title=element_text(size=9))+ggtitle(colnames(scores_2)[i])

  lay <- rbind(c(1,2,3),
               c(4,4,4))
  gridExtra::grid.arrange(p1,p2,p3,p4, layout_matrix = lay)
}
dev.off()



df_umap$a_gene <- GetAssayData()
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_gene))+geom_point(alpha=0.5,size=1)+
  scale_color_viridis(direction = -1)+theme_classic()+ #paste0("Varimax-", pc_num)
  theme(text = element_text(size=9), legend.title=element_text(size=9))+ggtitle(colnames(scores_2)[i])

rot_df$a_gene <-  scores_3[[i]]
ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=a_geneset))+geom_point(alpha=0.7,size=1.3)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_viridis(direction = -1, option='inferno')+theme_classic()+
  theme(text = element_text(size=9),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=9,angle=90,color='black'),
        legend.title = element_blank(),legend.position="top")+
  guides(colour=guide_colourbar(barwidth=15))+
  ggtitle(colnames(scores_2)[i])





pdf('TNFa_genesetUcell_varimax.pdf')
for(i in 1:5){
  rot_df$a_geneset <-  scores_3[[i]]
  p=ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=a_geneset))+geom_point(alpha=0.7,size=1.3)+
    theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
    scale_color_viridis(direction = -1, option='inferno')+theme_classic()+
    theme(text = element_text(size=18),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=13,angle=90,color='black'),
          legend.title = element_blank(),legend.position="top")+
    guides(colour=guide_colourbar(barwidth=15))+
    ggtitle(colnames(scores_2)[i])
  print(p)
  
  
  df_umap$a_geneset <-  scores_3[[i]]
  p=ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=a_geneset))+geom_point(alpha=0.5,size=1)+
    scale_color_viridis(direction = -1, option = 'inferno',name='Ucell TNF')+theme_classic()+ #paste0("Varimax-", pc_num)
    theme(text = element_text(size=22), legend.title=element_text(size=17))+ggtitle(colnames(scores_2)[i])
  print(p)
}
dev.off()


rownames(merged_samples)[grep('Tnf',rownames(merged_samples))]
rownames(merged_samples)[grep('Tnfsf1a',rownames(merged_samples))]



