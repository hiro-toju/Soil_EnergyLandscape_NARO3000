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
load.lib(c('ggplot2', 'tidyr'))

# -- Create directory to save
dir <- make.dir('04_Functions/04_output')

# -- Load data table
sml <- readRDS("Table/sample_info_01.rds") ; dim(sml)


picrust <- read.table('04_Functions/Table/picrust2Output2/pathways_out/path_abun_unstrat_descrip.tsv.gz', 
                      sep='\t', quote='', header=TRUE)
picrustT <- t(picrust)

fngtrait <- read.csv('04_Functions/Table/FungalTraits 1.2_ver_16Dec_2020 - V.1.2.csv')

## ========================================================== ##
