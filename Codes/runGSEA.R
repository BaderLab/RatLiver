### run GSEA analysis on the loadings.rnk files

# Ruth GSEA tutorial
# gsea-cli.sh GSEA [parameters]
# https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/Protocol2_createEM.html

# GSEA file format descriptions
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29


## in GSEA 4 the following parameters have been specified
# rnk - path to the rank file
# gmx - path to the gene set definition (gmt) file
# collapse - true/false indicates whether the expression/rnk file needs to be collapsed from probes to gene symbols
# nperm - number of permutations
# scoring_scheme -
#   rpt_label - name of the directory with output
# rnd_seed - random seed to use
# set_max - maximum size for individual gene sets. In GSEA interface this is set to 500 but we prefer to use a more stringent setting of 200.
# set_min - minimum size for individual gene sets
# zip_report - true/false to zip output directory
# out - directory where to place the result directory.


source('Codes/Functions.R')
Initialize()
gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'



####### Running GSEA on the markers list ########
rnk_file_path = paste0('~/RatLiver/Results/strain_variation_loadings/ranked_files/')
rnk_files <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = T)
working_dir = paste0('~/RatLiver/Results/strain_variation_loadings/gsea_results/')
dir.create(working_dir)

for(i in 1:length(rnk_files)){
  rnk_file = rnk_files[i]
  print(rnk_file)
  analysis_name = gsub('_loadings.rnk', '',list.files(rnk_file_path)[i])
  working_dir = working_dir

  GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                       rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                       analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , #set min=15
                       working_dir, " > gsea_output.txt")
  system(GSEA_command)
}




####### Running GSEA on the loadings list ########

source('Codes/Functions.R')
Initialize()
gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file = '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'

rnk_file_path = '~/RatLiver/Results/rot_PCs/'
rnk_file_path = '~/RatLiver/Results/new_samples/rotated_loadings/'
rnk_file_path = '~/RatLiver/Results/old_samples/rotated_loadings/'

#### selecting the subset of the files to run gsva on ####

#pc_indices <- c(5, 2, 15, 16, 3, 25, 20, 24, 11, 6, 4, 10, 8, 14, 9, 29, 12, 7)
map = 'set_1' #'set_2'
pc_indices <- c(9, 17, 29) #set1
if(map=='set_2') pc_indices <- c(9, 13) # set-2, last 3 need to be made


rnk_files_to_include <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = F) %in% 
  paste0('rot_PC',pc_indices,'_loadings.rnk')

#rnk_files <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = T)
rnk_files <- list.files(rnk_file_path,pattern = '*.rnk' ,full.names = T)[rnk_files_to_include]
######

working_dir = paste0('~/RatLiver/Results/new_samples/gsea_rotated_loadings') #gsea_rotated_loadings_immune
working_dir = paste0('~/RatLiver/Results/old_samples/gsea_rotated_loadings') 
dir.create(working_dir)

for(i in 1:length(rnk_files)){
  rnk_file = rnk_files[i]
  print(rnk_file)
  # analysis_name = gsub('.rnk', '',list.files(rnk_file_path)[rnk_files_to_include][i])
  analysis_name = gsub('.rnk', '',list.files(rnk_file_path)[rnk_files_to_include][i])
  print(analysis_name)
  working_dir = working_dir
  GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                       rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                       analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , #set min=15
                       working_dir, " > gsea_output.txt")
  system(GSEA_command)
}






