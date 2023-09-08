############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("/Users/toju/Dropbox/NARO_3000soil/Statistics")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))

## -- Create directory to save
dir <- make.dir('05_Community_structure/01_output')

# -- Load data table
sml <- readRDS("Table/0301_sml.rds") 

comm16s <- readRDS('Table/rrarefy_16S.rds')
commits <- readRDS('Table/rrarefy_ITS.rds')
dflist <- list(comm16s, commits)
dflist <- lapply(dflist, function(x){ x <- as.data.frame(x)[rownames(sml),]; rownames(x) <- rownames(sml); return(x)})

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxa16S[,-c(1,2)] <-  apply(taxa16S[,-c(1,2)], 2, function(x){paste(substr(taxa16S[,1], 1, 3),x,sep='_')})

taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxaITS[,-c(1,2)] <-  apply(taxaITS[,-c(1,2)], 2, function(x){paste('Fng',x,sep='_')})

taxa <- rbind(taxa16S, taxaITS)

############################################################################
# ASV

commonrow <- intersect(rownames(dflist[[1]]),rownames(dflist[[2]]))
bf <- na.omit(cbind(dflist[[1]][commonrow, ], dflist[[2]][commonrow, ]))
bf[bf>0] <- 1

reldf <- lapply(dflist, function(x){ na.omit(x/rowSums(x)) })
distdf <- lapply(reldf, vegdist, method='bray')
pcalist <- lapply(distdf, cmdscale)
#mdslist <- lapply(reldf, metaMDS, distance='bray', parallel=parallel::detectCores())
#tsnelist <- lapply(reldf, Rtsne, check_duplicates=FALSE)
#phate <- lapply(reldf, phate)

saveRDS(pcalist, sprintf('%s/dimensionReduction_rrarefy.rds', dir$rdsdir))

#saveRDS(list(pcalist,  tsnelist, phate), sprintf('%s/dimensionReduction_rrarefy.rds', dir$rdsdir))

#pcoa <- cmdscale(vegdist(bf, method='jaccard', binary=TRUE))
#mdslist <- metaMDS(bf, distance='jaccard', parallel=parallel::detectCores())
#tsnelist <- Rtsne(bf, check_duplicates=FALSE)
#phate <- phate(bf)
#saveRDS(list(pcoa, tsnelist, phate), sprintf('%s/dimensionReduction_binary_rrarefy.rds', dir$rdsdir))

############################################################################
