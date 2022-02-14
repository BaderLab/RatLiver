source('Codes/Functions.R')
Initialize()

################################
##### merging all the information related to the set-1
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


top_n = 30
get_head_genes <- function(df){df$genes[1:top_n] }
get_tail_genes <- function(df){ df$genes[(nrow(df)-top_n):nrow(df)] }

get_head_10genes <- function(df){df$genes[1:10] }
get_tail_10genes <- function(df){ df$genes[(nrow(df)-10):nrow(df)] }


old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
set1_data <- your_scRNAseq_data_object
set1_data@meta.data$cluster = as.character(sCVdata_list$res.0.6@Clusters)
set1_data@meta.data$strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2))

varimax_res_set1 <- readRDS('Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')
score_set1 <- data.frame(varimax_res_set1$rotScores)
colnames(score_set1) = paste0('Varimax_', 1:ncol(score_set1))
loading_set1 <- varimax_res_set1$rotLoadings

set1_info <- read.csv('figure_panel/set-1-final-info.csv')
colnames(set1_info)[1] = 'clusters'
set1_data@meta.data$umi = colnames(set1_data)
meta.data.label = merge(set1_data@meta.data, set1_info, by.x='cluster', by.y='clusters')
meta.data.label <- meta.data.label[match(colnames(set1_data),meta.data.label$umi),]
head(meta.data.label)
set1_data@meta.data = meta.data.label
table(set1_data$label)
set1_data$sample = unlist(lapply(str_split(set1_data$umi, '_'), function(x) paste0(x[2:3], collapse = '_')))

set1_df <- data.frame(umi=rownames(score_set1), 
                      strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2)), 
                      cluster = set1_data$cluster, 
                      meta.data.label,
                      score_set1)
head(set1_df)


############################ 
##### Frequency of lyz2 and Cd68 markers over CD45+ cells for each strain - check both set-1 and set2 maps
set1_data@meta.data$strain = unlist(lapply(str_split(colnames(set1_data),'_'), '[[', 2))

data = set1_data
marker.df = data.frame(strain= data$strain,
                       Ptprc = ifelse(GetAssayData(data)[rownames(data) == 'Ptprc',] > 0, 'Ptprc+', 'Ptprc-'),
                       Cd68 = ifelse(GetAssayData(data)[rownames(data) == 'Cd68',] > 0,'Cd68+', 'Cd68-'),
                       Lyz2 = ifelse(GetAssayData(data)[rownames(data) == 'Lyz2',] > 0,'Lyz2+', 'Lyz2-'),
                       Itgal = ifelse(GetAssayData(data)[rownames(data) == 'Itgal',] > 0,'Itgal+', 'Itgal-'),
                       Marco = ifelse(GetAssayData(data)[rownames(data) == 'Marco',] > 0,'Macro+', 'Macro-'))


marker.df$Markers = paste0(marker.df$Ptprc, marker.df$Cd68, marker.df$Itgal)
marker.df2 = data.frame(table(Strain=marker.df$strain, Markers=marker.df$Markers))
head(marker.df2)
ggplot(data=marker.df2, aes(x=Markers, y=Freq, fill=Strain)) +
  geom_bar(color='black', stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("pink1", "skyblue1"))+ylab('Number of cell in the Set-1 map')+
  theme_classic()+xlab('')+theme(axis.text.x = element_text(angle = 90, size=12))



#### Looking at Ptprc+ cells only 
marker.df3 = marker.df[marker.df$Ptprc=='Ptprc+',]

### evaluating double marker expression
marker.df3$Markers = paste0( marker.df3$Cd68, marker.df3$Itgal)
marker.df4 = data.frame(table(Strain=marker.df3$strain, Markers=marker.df3$Markers))

### evaluating single marker expression
marker.df4 = data.frame(table(Strain=marker.df3$strain, Markers=marker.df3$Itgal))

ggplot(data=marker.df4, aes(x=Markers, y=Freq, fill=Strain)) +
  geom_bar(color='black',stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  ylab('Number of Ptprc+ cell in the Set-1 map')+
  theme_classic()+xlab('')+
  theme(axis.text.x = element_text(angle = 90, size=12))


###################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


marker.df_bin = marker.df
marker.df_bin$sample = unlist(lapply(str_split(data$umi, '_'), function(x) paste0(x[2:3], collapse = '_')))

marker.df_bin$Ptprc.Cd68 = paste0(marker.df_bin$Ptprc, marker.df_bin$Cd68)
marker.df_bin$Ptprc.Marco = paste0(marker.df_bin$Ptprc, marker.df_bin$Marco)
marker.df_bin$Ptprc.Itgal = paste0(marker.df_bin$Ptprc, marker.df_bin$Itgal)

#marker.df_bin$Ptprc.Lyz2 = paste0(marker.df_bin$Ptprc, marker.df_bin$Lyz2)

head(marker.df_bin)
marker.df_bin.l=data.frame(markers=c(marker.df_bin$Ptprc.Cd68, 
                                     marker.df_bin$Ptprc.Marco, 
                                     marker.df_bin$Ptprc.Itgal),
                           strain = rep(marker.df_bin$strain,3),
                           sample = rep(marker.df_bin$sample,3))
table(marker.df_bin.l$markers)
marker.df_bin.l = marker.df_bin.l[marker.df_bin.l$markers %in% c('Ptprc+Cd68+', 'Ptprc+Itgal+', 'Ptprc+Macro+'),]
head(marker.df_bin.l)

tab = table(marker.df_bin.l$markers,
            marker.df_bin.l$sample, 
            marker.df_bin.l$strain)
tab = data.frame(rbind(melt(tab[,,1]),melt(tab[,,2])))
colnames(tab) = c('markers','sample', 'freq')
tab$strain = unlist(lapply(str_split(tab$sample, '_'), '[[', 1))
tab = tab[tab$freq != 0,]

ggplot(data=tab, aes(x=markers, y=freq,fill=sample)) +
  geom_bar(color='black',stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  ylab('Cell count')+
  theme_classic()+xlab('')

tab2 = data_summary(tab[,-2], 'freq', c('markers','strain'))
ggplot(data=tab2, aes(x=markers, y=freq,fill=strain)) +
  geom_bar(color='black',stat="identity", position=position_dodge())+
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  ylab('Cell count')+
  theme_classic()+xlab('')+
  geom_errorbar(aes(ymin=freq-sd, ymax=freq+sd), width=.2,
                position=position_dodge(.9)) 


#### set1 - significant relationship between strain and Cd68 in the Ptprc+ cells - p.val=1.215e-07
#### set1 - significant relationship between strain and Lyz2 in the Ptprc+ cells - p.val= 0.007681

#### set2 - significant relationship between strain and Cd68 in the Ptprc+ cells - p.val=1.067e-12
#### set2 - NO significant relationship between strain and Lyz2 in the Ptprc+ cells - p.val= 0.7605

conting_table = table(marker.df3$strain, marker.df3$Cd68) 
conting_table = table(marker.df3$strain, marker.df3$Lyz2)

fisher.test(conting_table)
chisq.test(conting_table)


###### check these differrences in cluster number-5 which includes the non-Inf Mac 
# do you think these differences might be related to varriations in number of cells ? sample imbalance
# if donw-sample, should we start from raw data? prenormalization might effect the results??


################ checking histology statistics
CD3_data = read.csv('Data/Histology/CD3.csv')
CD68_data = read.csv('Data/Histology/CD68.csv')
CD8_data = read.csv('Data/Histology/CD8.csv')

his_data = CD8_data
rownames(his_data) = paste0('layer_',his_data[,1])
his_data = his_data[,-1]
dim(his_data)
his_data_DA = his_data[,1:30]
his_data_LEW = his_data[,31:ncol(his_data)]

#### total differeces between the DA and LEW  
t.test(colSums(his_data_DA), colSums(his_data_LEW))
wilcox.test(colSums(his_data_DA), colSums(his_data_LEW))

#### CD3 general p-value: 0.03442 is significant
#### CD68 general p-value: 0.2374 is Not significant
#### CD8 general p-value: 0.3881 is Not significant

#### Layer-specific differeces between the DA and LEW ### alternative test
t.test.df = data.frame(layer=rownames(his_data), p.value=rep(0,10))
wilcoxon.test.df = data.frame(layer=rownames(his_data), p.value=rep(0,10))

for(i in 1:nrow(his_data_DA)){
  t.test.df[i,2] = t.test(his_data_DA[i,], his_data_LEW[i,])$p.value
  wilcoxon.test.df[i,2] = t.test(his_data_DA[i,], his_data_LEW[i,])$p.value
}
t.test.df$sig = t.test.df$p.value < 0.05/10
t.test.df$sig

wilcoxon.test.df$sig = wilcoxon.test.df$p.value < 0.05/10
wilcoxon.test.df$sig
#### CD3 layer-9 p-value: 0.002447611 is significant
#### CD68 not sig
#### CD8 not sig
