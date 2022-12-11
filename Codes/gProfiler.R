source('Codes/Functions.R')
Initialize()


get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),
                  organism = model_animal_name)
  return(gostres)
}


#### importing the input data
INPUT_NAME = 'rat_Rnor'  # 'rat_Rnor''rat_Lew_01', 'rat_DA' 'mouse' 
model_animal_name = "rnorvegicus" #'mmusculus' 

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)

RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
if(INPUT_NAME=='rat_Lew_01')  RES = '_res.0.8' 

Cluster_markers_merged <- readRDS('~/Desktop/cell_type_markers/rat_Rnor/markers_rat_Rnor_mito_50_lib_1500_res.1.rds')
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME, '/markers/markers_', 
                                         OUTPUT_NAME, RES,'.rds'))
cluster_names <- names(Cluster_markers_merged)
lapply(Cluster_markers_merged, dim)

######## enrichment using gProfiler
gp_enrich_res <- sapply(1:length(Cluster_markers_merged), 
                        function(i){
                          markers = Cluster_markers_merged[[i]]$ensemble_ids
                          enrich_res = get_gprofiler_enrich(markers, model_animal_name);return(enrich_res)
                        },simplify = F)
names(gp_enrich_res) <- names(Cluster_markers_merged)

lapply(gp_enrich_res, function(x) head(x$result))
saveRDS(gp_enrich_res,'~/Desktop/cell_type_markers/rat_Rnor/gProfiler_rat_Rnor_mito_50_lib_1500_res.1.rds')


gp_enrich_res_filt <- lapply(gp_enrich_res, 
                             function(x) x$result[x$result$query_size>5 & x$result$query_size<350 & x$result$intersection_size >3, ])
names(gp_enrich_res_filt)  
lapply(gp_enrich_res_filt, head)
head(gp_enrich_res_filt[['cluster_10']], 10)

p <- gostplot(gostres, capped = T, interactive = FALSE)
# pt2 <- publish_gosttable(gostres, use_colors = TRUE, filename = NULL)
