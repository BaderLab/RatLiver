
library(tidyverse)
BiocManager::install('hrbrthemes')
library(hrbrthemes)
library(viridis)
library(reshape2)

mac_hist = data.frame(t(read.csv('CD68.csv')))
colnames(mac_hist) = paste0('L', 1:10)
mac_hist = mac_hist[-1,]
mac_hist$strain = unlist(lapply(strsplit(rownames(mac_hist), '\\.'), '[[', 1))
mac_hist.df = melt(mac_hist, id.vars = "strain")
colnames(mac_hist.df)[2] = 'Layer'
head(mac_hist.df)
ggplot(mac_hist.df, aes(fill=strain, y=value, x=Layer)) + 
  geom_bar(position="dodge", stat="identity")+theme_classic()+geom_jitter(color="black", size=0.4, alpha=0.9)

