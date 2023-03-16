source('Codes/Functions.R')
source('Codes/convert_human_to_ortholog_functions.R')
Initialize()
library(SoupX)
library(Matrix)

# note: DA_M_10WK_004_strained is the 200515_A00827_0155_AHLV55DRXX_McParland_Sonya reseq sample
sample_names = list.dirs('~/RatLiver/Data/SoupX_data/SoupX_inputs/',recursive = FALSE, full.names = FALSE)
i = 1


rho_list = sapply(1:length(sample_names), function(i){
  sample_name = sample_names[i]
  print(sample_name)
  soupX_out = readRDS(paste0('~/RatLiver/Data/SoupX_data/SoupX_outputs/default_param/', sample_name, '_soupX_out.rds'))
  print(table(soupX_out$sc$metaData$rho))
  return(soupX_out$sc$metaData$rho[i])
})

names(rho_list) = sample_names

add_val = 0.1
for(i in 1:length(sample_names)){
  sample_name = sample_names[i]
  print(sample_name)
  sc = load10X(paste0('~/RatLiver/Data/SoupX_data/SoupX_inputs/', sample_name, '/'))
  
  #sc = autoEstCont(sc, forceAccept = TRUE)
  sc = setContaminationFraction(sc,rho_list[i]+add_val, forceAccept = TRUE)
  
  out = adjustCounts(sc)
  saveRDS(list(sc=sc,out=out),paste0('~/RatLiver/Data/SoupX_data/SoupX_outputs/', sample_name, '_soupX_out_Rhoplus.',add_val,'.rds'))
}
dev.off()


sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
sc = estimateSoup(sc)

soupProf = data.frame(row.names = rownames(toc),
                      est = rowSums(toc)/sum(toc),
                      counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops,soupProf)


sc = setClusters(sc,setNames(PBMC_metaData$Cluster,rownames(PBMC_metaData)))





for(i in 1:length(sample_names)){
  sample_name = sample_names[i]
  print(sample_name)
  sc = load10X(paste0('~/RatLiver/Data/SoupX_data/SoupX_inputs/', sample_name, '/'))
  sc = autoEstCont(sc, forceAccept = TRUE)
  out = adjustCounts(sc)
  #saveRDS(list(sc=sc,out=out),paste0('~/RatLiver/Data/SoupX_data/SoupX_outputs/', sample_name, '_soupX_out.rds'))
}
