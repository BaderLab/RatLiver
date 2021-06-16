devtools::install_github("BaderLab/CCInx")
library(CCInx)
#> Loading required package: shiny
load(system.file("DemoData/DemoDE.RData",package="CCInx"))
lapply(deL,head)

## Ranking nodes by differential expression between conditions

inx <- BuildCCInx(GeneStatList=deL,
                  GeneMagnitude="logFC",
                  GeneStatistic="padj",
                  Species="mmusculus")
head(inx$edges)
head(inx$nodes)
PlotCCInx(INX=inx,
          cellTypeA="GABA",cellTypeB="DOPA",
          proteinTypeA="Receptor",proteinTypeB="Ligand",
          TopEdges=50)
# ViewCCInx(inx) ## shiny viewer
load(system.file("DemoData/DemoExpr.RData",package="CCInx"))
show(e13cortex)
gsl <- BuildGeneStatList(inD=e13cortex,
                         cl=colData(e13cortex)$cellTypes,
                         assayType="logcounts")
lapply(gsl[1:3],head)

inx <- BuildCCInx(GeneStatList=gsl,
                  Species="mmusculus")

head(inx$edges)
head(inx$nodes)
PlotCCInx(INX=inx,
          cellTypeA="ProjectionNeurons",cellTypeB="CorticalPrecursors",
          proteinTypeA="Ligand",proteinTypeB="Receptor",
          GeneMagnitudeThreshold=.5)




