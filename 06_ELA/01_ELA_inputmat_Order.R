############################################################################
####
#### R script for Fujita (2022)
####
#### Energy landscape analysis
#### 2022. 3. 16 Fujita
#### R 4.1.2
#### setwd("~/Desktop/NARO/220520/")
#### 
############################################################################

set.seed(123)
setwd("../")

lv <- "Order"

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('ggplot2', 'tidyr','cowplot','Rtsne','vegan','RColorBrewer','phateR'))
library(mgcv)
library(rgl)

# -- Create directory to save
dir <- make.dir('06_ELA/01_output')

############################################################################
# -- Load community matrix data

rarepro <- readRDS('Table/rrarefy_16S.rds')
rarefun <- readRDS('Table/rrarefy_ITS.rds')

bothrow <- intersect(rownames(rarepro), rownames(rarefun))

############################################################################
# Selecting data

sml <- readRDS("Table/0301_sml.rds") 

#sml2 <- sml[, c('sampling_date', 'latitude', 'longitude', 'crop', 'former_crop', 'DL', 'pH_dry_soil', 'EC_electric_conductivity', 'CNratio', 'available_P')]

env <- na.omit(sml[, c('pH_dry_soil', 'EC_electric_conductivity', 'CNratio', 'available_P')])

rl <- intersect(bothrow, rownames(env))

saveRDS(sml[rl, ], sprintf('Table/sml_%s.rds', length(rl)))
saveRDS(env[rl, ], sprintf('Table/env_%s.rds', length(rl))) # to be standardized before ELA.

############################################################################

inpro <- rarepro[rl, ]
infun <-  rarefun[rl,]

inpro <- inpro[, colSums(inpro)>0]
infun <- infun[, colSums(infun)>0]
inboth <- cbind(inpro, infun)

dim(inpro)
dim(infun)
dim(inboth)

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxa <- rbind(taxa16S, taxaITS)

############################################################################
# 
n.core <- parallel::detectCores()

txpro <-  Taxa.mat(inpro, taxa=taxa, taxaLabel=lv)
txfun <-  Taxa.mat(infun, taxa=taxa, taxaLabel=lv)
txboth <-  Taxa.mat(inboth, taxa=taxa, taxaLabel=lv) # "Unidentified" of pro and fun data are merged. However, they will be removed in the taxa screening process before ELA. 

saveRDS(txpro, sprintf('Table/ELA_matrix_pro_%s.rds', lv))
saveRDS(txfun, sprintf('Table/ELA_matrix_fun_%s.rds', lv))
saveRDS(txboth, sprintf('Table/ELA_matrix_both_%s.rds', lv))

dim(txpro)
dim(txfun)
dim(txboth)

############################################################################

