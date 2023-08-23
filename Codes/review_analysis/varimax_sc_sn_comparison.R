varimax_snuc = readRDS('~/rat_sham_sn_data/standardQC_results/sham_sn_merged_standardQC_varimax_res2.5_updated_July23.rds')
varimax_sc <- readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds') ## MT-removed

rotatedLoadings.snuc <-data.frame(sapply(1:ncol(varimax_snuc$loading), function(i)varimax_snuc$loading[,i]))
colnames(rotatedLoadings.snuc) = paste0('Var ',c(1:ncol(rotatedLoadings.snuc)), '   (sn)')
rotatedLoadings.snuc <- rotatedLoadings.snuc[,1:16]
rotatedLoadings.snuc$gene = rownames(rotatedLoadings.snuc)

rotatedLoadings.sc <-data.frame(sapply(1:ncol(varimax_sc$rotLoadings), function(i)varimax_sc$rotLoadings[,i]))
colnames(rotatedLoadings.sc) = paste0('Var ',c(1:ncol(rotatedLoadings.sc)), '   (sc)')
rotatedLoadings.sc <- rotatedLoadings.sc[,1:15]
rotatedLoadings.sc$gene = rownames(rotatedLoadings.sc)

merged_var_df = merge(rotatedLoadings.sc, rotatedLoadings.snuc, by.x='gene', by.y='gene')
dim(merged_var_df)
col.sn = colnames(merged_var_df)[grep(colnames(merged_var_df), pattern = 'sn')]
col.sc = colnames(merged_var_df)[grep(colnames(merged_var_df), pattern = 'sc')]

cor_varimax = cor(merged_var_df[,-1])
cor_varimax.sub = cor_varimax[colnames(cor_varimax) %in% col.sn, row.names(cor_varimax) %in%col.sc]
dim(cor_varimax.sub)
pheatmap(cor_varimax.sub, fontsize = 12)



var16.df = data.frame(Gene=rotatedLoadings.snuc$gene,score=rotatedLoadings.snuc$`Var 16   (sn)`)  
var16.df = var16.df[order(var16.df$score, decreasing = T),]
row.names(var16.df) <- NULL
var16.df$score   = round(var16.df$score ,3)
var16.df <- var16.df[!grepl(var16.df$Gene, pattern = 'Mt'),]
colnames(var16.df)[2] = 'Loading score'
dev.off()
gridExtra::grid.table(head(var16.df,10))
head(var16.df,20)


embedd_df_rotated = data.frame(varimax_snuc$score)
pc_num = 16
rot_df <- data.frame(Varimax_1=embedd_df_rotated$Var.PC_1,
                     emb_val=embedd_df_rotated[,pc_num+1],
                     label=as.character(embedd_df_rotated$labels),
                     #Library_size=merged_samples$nCount_RNA,
                     #num_expressed_genes=merged_samples$nFeature_RNA,
                     strain=embedd_df_rotated$strain,
                     #label=factor(df_umap$label, levels = as.character(set_1_cl_ord)),
                     sample_name = embedd_df_rotated$sample_name
)

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=label))+geom_point(alpha=0.7,size=1)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(name='Cell Type',values = Colors)+theme_classic()+
  theme(text = element_text(size=18))#, legend.title = element_blank()

  #ggtitle(paste0('Set-2 Immune-subclusters over Varimax 1 and ', pc_num))
ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.7,size=3)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  theme(text = element_text(size=18))#, legend.title = element_blank()

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point(alpha=0.7,size=3)+
  theme_classic()+ylab(paste0('Varimax ',pc_num))+xlab("Varimax 1")+
  scale_color_brewer(palette = "Set3")+theme_classic()+
  theme(text = element_text(size=18))#, legend.title = element_blank()


#### boxplots plots for varimax
ggplot(rot_df, aes(x=label, y=emb_val, fill=label))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = mycolors)+theme_classic()+å
  +ylab(paste0('Varimax ', pc_num))+
  theme(text = element_text(size=18),
        axis.text.x = element_text(size=13,angle=90,color='black'),
        legend.title = element_blank(), legend.position = "none")+
  #scale_x_discrete(labels = xlabel2)+
  xlab('')

ggplot(rot_df, aes(x=label, y=emb_val, fill=strain))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  ylab(paste0('Varimax ', pc_num))+
  theme(text = element_text(size=18),
        axis.text.x = element_text(size=13,angle=90,color='black'),
        legend.title = element_blank(), legend.position = "none")+
  scale_x_discrete(labels = xlabel2)+xlab('')


df_umap$Varimax = abs(rot_df$emb_val)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=abs(Varimax)))+geom_point(alpha=0.5,size=1)+
  scale_color_viridis(direction = -1, option = 'inferno',name=paste0("Varimax ", pc_num))+theme_classic()+
  theme(text = element_text(size=22), legend.title=element_text(size=17))
#+ggtitle(paste0('Varimax ', pc_num, ' scores over set2\nimmune-subclusters UMAP'))
#'plasma', 'inferno', 'magma'

###### strain-specific varimax factors

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=strain))+geom_point(alpha=0.5,size=1.5)+
  theme_classic()+ylab(paste0('varimax_',pc_num))+xlab("varimax 1")+
  scale_color_manual(values = c("#CC79A7","#0072B2"))+theme_classic()+
  theme(text = element_text(size=16), legend.title = element_blank())# "#56B4E9"
#+ggtitle(paste0('Distribution of cells based on strain\nover Varimax ', pc_num))

###### strain-specific varimax factors - sample

ggplot(rot_df, aes(x=Varimax_1, y=emb_val, color=sample_name))+geom_point(alpha=0.5,size=1.4)+
  theme_classic()+ylab(paste0('Varimax-',pc_num))+xlab("Varimax-1")+theme_classic()+
  theme(text = element_text(size=18), legend.title = element_blank())# "#56B4E9"
#+ggtitle(paste0('Distribution of cells based on strain\nover Varimax ', pc_num))
