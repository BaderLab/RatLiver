

merged_samples <- readRDS('cell_browser/snRNAseq_ratLiver_cellBrowser.rds')
head(merged_samples@meta.data )
table(merged_samples$annotation)
ncol(merged_samples)
colnames(merged_samples)

Idents(merged_samples) = (merged_samples$annotation)
merged_samples_sub = merged_samples[,!merged_samples$annotation %in% c('Hep 0', 'Hep 1','Hep 2','Hep 3','Unknown/High Mito')]
merged_samples
table(merged_samples_sub$annotation)
Cholagiocytes = 29
Cholagiocytes_Markers = c('Epcam', 'Sox9', 'Anxa4', 'Krt8')


General_Macrophage_Markers = c('Cd68', 'Clec4f')

Non_inf_mac = 19 
Non_inf_mac_Markers = c('Cd5l', 'Marco', 'Cd163', 'C1qa', 'C1qc', 'Aif1', 
                        'Hmox1', 'Vsig4',  'Ccl6') #'Slc11a' Ptas1

Inf_mac = 33 
Inf_mac_Markers = c('Cd74', 'RT1-Ba', 'RT1-Bb', 'RT1-Da', 'RT1-Db1', 'Lyz2', 'Ighm')


Mesenchymal = c(24) 
Mesenchymal_Markers = c('Col3a1', 'Colec10', 'Colec11', 'Ecm1', 
                        'Calcrl' , 'Col1a1', 'Col6a2' , 'Lrat', 'Reln', 'Hgf', 'Pth1r', 'Tagln') #'Stab2', 'Vwf', 


Endothelial = c(11, 30) 
Endothelial_Markers = c('Aqp1', 'Fcgr2b', 'Gpr182', 'Lyve1', 'Stab2', 'Bmp2', 'Ramp2', 'Ctsl', 'Stab2', 'Sparc', 'Eng' )


clusters_to_include = c(CV_Hep, other_Hep, Periportal_Hep, Cholagiocytes, 
                        Non_inf_mac, Inf_mac, Mesenchymal, Endothelial)

markers_list = c( Cholagiocytes_Markers, General_Macrophage_Markers, Non_inf_mac_Markers, Inf_mac_Markers, 
                 Mesenchymal_Markers, Endothelial_Markers) 

Idents(merged_samples_sub) = merged_samples_sub$cluster
table(merged_samples_sub$cluster)
markers_list = unique(markers_list)

nucseq_ord = c("Cholangiocytes (29)", "Non-inf Macs (19)" ,"Inf Macs (33)",
"Mesenchymal (24)"  , "Endothelial (11)", "Endothelial (30)")

Idents(merged_samples_sub) <- factor(merged_samples_sub$cluster, levels = nucseq_ord)
DotPlot(merged_samples_sub, features = markers_list) + RotatedAxis()+ xlab('Markers')+
  ylab('')+theme(axis.text.x = element_text(size=12))

############## Dot plots ###############
############### Set-1 ########
merged_samples2 = merged_samples
dim(merged_samples2)
set1_info_ord = set1_info[match(set_2_cl_ord, set1_info$label),]
matrix = as.matrix(set1_info_ord[,4:8])
markers = c()
for(i in 1:nrow(set1_info_ord)) markers = c(markers,matrix[i,])
markers <- markers[markers!='']
markers = unique(unname(markers))
div_num = 42#34 #58 #64
markers_final = unique(c(markers[1:div_num],'Ptprc', 'Cd68', markers[(div_num+1):length(markers)]))
markers_final = unique(c('Ptprc', 'Cd68',markers))
markers_final = unique(c('Ptprc', markers[1:div_num], 'Cd68', markers[(div_num+1):length(markers)]))



markers_final=c("Cps1","Scd","Cyp2a1","Itih3","Fgb","Calr","C3","Itih4","Arg1","G6pc","Slco1b2","Igfbp1",
                "Insig1","Hamp","Fabp1","Dbi","Apoa2","Cox8a","Alb","Apoc1","Fth1","Apoc3",  
                "Cox7b","Rpl10",'Ifi27',"Ctsl","Lyve1","Fcgr2b","Sparc",
                'Rspo3', 'Clec4g', 'Bmp2',
                "Stab2","Fam167b","Eng",'Ramp2',"Calcrl","Ecm1","Col3a1","Igfbp7",
                "Ptprc","Cd68","Clec4f",
                "Marco","Vsig4","Cd5l","Hmox1","Cd163","Cd74","Lyz2",   
                "RT1-Db1", "RT1-Da","Il18","Klrd1","Gzmk","Cd7","Ccl5","Gzma")



Idents(merged_samples2) <- factor(df_umap$label, levels = set_2_cl_ord)
DotPlot(merged_samples2, features = markers_final) + RotatedAxis()+ xlab('Markers')+
  ylab('')+theme(axis.text.x = element_text(size=12))
#split.by = "groups"
DotPlot(merged_samples2, features = markers_final) + RotatedAxis()+ xlab('Markers')+
  ylab('')+theme(axis.text.x = element_text(size=12))
#split.by = "groups"
names(table(df_umap$label))
