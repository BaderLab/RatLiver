mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
## view the first few rows of the data
head(mydata)
summary(mydata)

xtabs(~admit + rank, data = mydata)
mydata$rank <- factor(mydata$rank)
mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")
summary(mylogit)


old_data_scClustViz_object <- "Results/old_samples/for_scClustViz_mergedOldSamples_mt40_lib1500_MTremoved.RData"
load(old_data_scClustViz_object)
cluster.df <- data.frame(cluster=sCVdata_list$res.0.6@Clusters)
cluster.df$cluster <- paste0('cluster_', as.character(cluster.df$cluster))
aFeature = 'cluster_3'
cluster.df$feature <- ifelse(cluster.df$cluster == aFeature, aFeature, 'other')
head(cluster.df)

varimax_df <- data.frame(readRDS('~/RatLiver/Results/old_samples/varimax_rotated_OldMergedSamples_mt40_lib1500_MTremoved.rds')$rotScores)
colnames(varimax_df) <- paste0('varPC_', 1:ncol(varimax_df))
varimax_df$cluster <- as.factor(cluster.df$feature)
varimax_df$strain <-as.factor(sapply(strsplit(x = rownames(varimax_df), '_'), '[[', 2))
table(varimax_df$strain)



## logistic regression does not work on data sets in which a feature perfectly seperates the datasets  
head(varimax_df)
varimax_df <- varimax_df[,-5]
mylogit <- glm(strain ~ . , data = varimax_df, family = "binomial", control=glm.control(maxit=50)) #+ varPC_2 + varPC_3
summary(mylogit)


### running random forrest on the data
##https://www.r-bloggers.com/2018/01/how-to-implement-random-forests-in-r/

#install.packages("randomForest")
library(randomForest)

varimax_df2 = varimax_df[,c(1:25,ncol(varimax_df))]
head(varimax_df2)
model0 <- randomForest(cluster ~ . , data = varimax_df2, importance = TRUE)
model1 <- randomForest(strain ~ . , data = varimax_df2, importance = TRUE)

importance(model1)        
importance(model0)        
varImpPlot(model1)        
varImpPlot(model0)        


