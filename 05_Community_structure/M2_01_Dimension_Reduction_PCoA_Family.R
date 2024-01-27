############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### 2023.12.21 Toju
#### R 4.3.2
#### 
############################################################################

set.seed(123)
setwd("../")

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))
library(foreach)
library(doParallel)
library(parallel)

#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)

lv <- "Family"

## -- Create directory to save
dir <- make.dir('05_Community_structure/M2_01_output')

# -- Load data table
rarepro <- readRDS('Table/rrarefy_16S.rds')
rarefun <- readRDS('Table/rrarefy_ITS.rds')

bothrow <- intersect(rownames(rarepro), rownames(rarefun))

############################################################################
# Selecting data

sml <- readRDS("Table/0301_sml.rds") 
env <- na.omit(sml[, c('pH_dry_soil', 'EC_electric_conductivity', 'CNratio', 'available_P')])

pl <- intersect(rownames(rarepro), rownames(env))
fl <- intersect(rownames(rarefun), rownames(env))

saveRDS(sml[pl, ], sprintf('Table/sml_pro_%s.rds', length(pl)))
saveRDS(env[pl, ], sprintf('Table/env_pro_%s.rds', length(pl))) # to be standardized before ELA.

saveRDS(sml[fl, ], sprintf('Table/sml_fun_%s.rds', length(fl)))
saveRDS(env[fl, ], sprintf('Table/env_fun_%s.rds', length(fl))) # to be standardized before ELA.

############################################################################

inpro <- rarepro[pl, ]
infun <-  rarefun[fl,]

inpro <- inpro[, colSums(inpro)>0]
infun <- infun[, colSums(infun)>0]

############################################################################
# Family

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')
taxa <- rbind(taxa16S, taxaITS)

txpro <-  Taxa.mat(inpro, taxa=taxa, taxaLabel=lv)
txfun <-  Taxa.mat(infun, taxa=taxa, taxaLabel=lv)

dflist <- list(txpro, txfun)

############################################################################
#pcalist <- lapply(dflist, prcomp, scale.=TRUE, rank.=2)

reldf <- lapply(dflist, function(x){ na.omit(x/rowSums(x)) })
distdf <- lapply(reldf, vegdist, method='bray')
pcalist <- lapply(distdf, cmdscale, eig=TRUE)
#pcalist <- lapply(distdf, cmdscale, k=2, eig=TRUE)

saveRDS(pcalist, sprintf('%s/PCoA_Family_rrarefy.rds', dir$rdsdir))
saveRDS(pcalist, 'Table/PCoA_Family_rrarefy.rds')


sink(sprintf('%s/PCoA_pro.txt', dir$tabledir))
eig_perc(pcalist[[1]]$eig, positive = T, plot = T)
dim(inpro)
sink()

sink(sprintf('%s/PCoA_fun.txt', dir$tabledir))
eig_perc(pcalist[[2]]$eig, positive = T, plot = T)
dim(infun)
sink()

############################################################################
