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

#bothrow <- intersect(rownames(rarepro), rownames(rarefun))

############################################################################
# Selecting data

sml <- readRDS("Table/0301_sml.rds") 

#sml2 <- sml[, c('sampling_date', 'latitude', 'longitude', 'crop', 'former_crop', 'DL', 'pH_dry_soil', 'EC_electric_conductivity', 'CNratio', 'available_P')]

env <- na.omit(sml[, c('pH_dry_soil', 'EC_electric_conductivity', 'CNratio', 'available_P')])

pl <- intersect(rownames(rarepro), rownames(env))
fl <- intersect( rownames(rarefun), rownames(env))

saveRDS(sml[pl, ], sprintf('Table/sml_pro_%s.rds', length(pl)))
saveRDS(env[pl, ], sprintf('Table/env_pro_%s.rds', length(pl))) # to be standardized before ELA.

saveRDS(sml[fl, ], sprintf('Table/sml_fun_%s.rds', length(fl)))
saveRDS(env[fl, ], sprintf('Table/env_fun_%s.rds', length(fl))) # to be standardized before ELA.

############################################################################

taxa16S <- readRDS('Table/taxonomyList_Prok.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxa <- rbind(taxa16S, taxaITS)

############################################################################
# Prokaryotes

inpro <- rarepro[pl, ]
inpro <- inpro[, colSums(inpro)>0]

n.core <- parallel::detectCores()

txpro <-  Taxa.mat(inpro, taxa=taxa, taxaLabel=lv)
saveRDS(txpro, sprintf('Table/ELA_matrix_each_pro_%s.rds', lv))

dim(txpro)

############################################################################
# Fungi

infun <- rarefun[fl, ]
infun <- infun[, colSums(infun)>0]

n.core <- parallel::detectCores()

txfun <-  Taxa.mat(infun, taxa=taxa, taxaLabel=lv)

saveRDS(txfun, sprintf('Table/ELA_matrix_each_fun_%s.rds', lv))

dim(txfun)
