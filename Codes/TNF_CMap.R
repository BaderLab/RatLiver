### Gary was referring to was Connectivity Map (https://clue.io/)
## which genes are transcribed in response to perturbations such as TNF treatment.
## full TNF response data from CMap
# matrix of normalized Z-scores from CMap for samples treated with TNF.  
# Column names include cell line, duration of treatment, and dosage for each sample.  Row names are gene symbols.  
# I've sorted it so highest average response are in the top rows.


load('~/RatLiver/TNFresponse.RData')
head(TNFlvl5)
average_cellLines = rowSums(TNFlvl5)/ncol(TNFlvl5)
hist(average_cellLines)
threshold = 1.5
sum(average_cellLines>threshold)
petrurb_tns = names(average_cellLines)[average_cellLines>threshold]
petrurb_tns.rat = .getMapped_hs2model_df(ensembl=ensembl, 
                                     candidateGenes=petrurb_tns, 
                                     model_animal_name=model_animal_name)

write.csv(petrurb_tns.rat, 'tnf_CMap_purturb_genes.csv')
petrurb_tns.rat = read.csv('tnf_CMap_purturb_genes.csv')
getUnemptyList(petrurb_tns.rat$rnorvegicus_homolog_associated_gene_name)


###### using the TNF secretion genes used in the paper -> doi: 10.3389/fimmu.2014.00538
model_animal_name = 'rnorvegicus' # 'mmusculus'

listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

lps_tnf_df = read.csv('LPS_TNF_geneset.csv')
lps_tnf.rat = .getMapped_hs2model_df(ensembl=ensembl, 
                       candidateGenes=lps_tnf_df$human, 
                       model_animal_name=model_animal_name)

lps_tnf_df = merge(lps_tnf_df, lps_tnf.rat, by.x='human', by.y='symbol', all.x=T)
write.csv(lps_tnf_df, 'lps_tnf_genes.csv')

##################
lps_tnf_df = read.csv('lps_tnf_genes.csv')
tnf_secretion.geneset = getUnemptyList(lps_tnf_df$rnorvegicus_homolog_associated_gene_name[lps_tnf_df$TNF=='+'])




