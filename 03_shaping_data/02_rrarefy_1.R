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

th <- 1000

## -- Loading Function and Library
library('AnalysisHelper')
load.lib(c('vegan'))

# -- Create directory to save
dir <- make.dir('03_shaping_data/02_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table
read16s <- as.data.frame(readRDS('Table/read_Prok.rds'))
readits <- as.data.frame(readRDS('Table/read_Fng.rds')) 
readlist <- list(Prok=read16s, Fng=readits)

sml <- readRDS("Table/0301_sml.rds") ; dim(sml)

################################################################
## rrarefy 

rrarefied16s <- rrarefy(read16s, th)
filt16s <- t(subset(t(rrarefied16s), rowSums(t(rrarefied16s)) > 0))
mean(1-rareslope(read16s, 1000))
pdf(sprintf("%s/rarecurve_16S.pdf", dir$figdir), w=5,h=5)
rarecurve(filt16s[sample(1:nrow(filt16s), 100),]); dev.off()

saveRDS(filt16s, file= sprintf("%s/rrarefy_16S.rds", dir$rdsdir))
Sample.ID <- rownames(filt16s)
out16s <- cbind(Sample.ID, filt16s)
#write.table(out16s, file="./02_output/Table/rrarefy_16S.txt", sep="\t", quote=F, row.names=F)


rrarefiedits <- rrarefy(readits, th)
filtits <- t(subset(t(rrarefiedits), rowSums(t(rrarefiedits)) > 0))
mean(1-rareslope(readits, 1000))
pdf(sprintf("%s/rarecurve_ITS.pdf", dir$figdir), w=5,h=5)
rarecurve(filtits[sample(1:nrow(filtits), 100),]); dev.off()

saveRDS(filtits, file=sprintf("%s/rrarefy_ITS.rds", dir$rdsdir))
Sample.ID <- rownames(filtits)
outits <- cbind(Sample.ID, filtits)
#write.table(outits, file="./02_output/Table/rrarefy_ITS.txt", sep="\t", quote=F, row.names=F)

dim(filt16s)
dim(filtits)
