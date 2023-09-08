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

# -- Create directory to save
library(AnalysisHelper)
library(seqinr)
dir <- make.dir('04_Functions/01_output')

# -- Load data table

dflist <- readRDS('Table/0303_rarefiedMatrix_rrarefy.rds')
taxaITS <- readRDS('Table/taxonomyList_Fng.rds')

taxa16S <- readRDS('Table/07_prokaryote_annotation.rds') 
seq16s <- read.fasta('Table/04_OTUseq_1.fasta') 
read16s <- readRDS('Table/07_prokaryote_seqtab.rds')[rownames(dflist[[1]]),colnames(dflist[[1]])]

############################################################################

## -- For picrust2
sub16s <- taxa16S[which(taxa16S[,1]%in% colnames(read16s)), ]
seqsub <- seq16s[rownames(sub16s)]
write.fasta(seqsub, sub16s[,1], sprintf("%s/picrutsInputFasta.fasta", dir$tabledir))

t16s <- t(read16s)
write.table(cbind(ASV=rownames(t16s), t16s),
            sprintf("%s/picrutsInput.txt", dir$tabledir), 
            row.names=FALSE, quote=FALSE, sep="\t")

## -- For FAPROTAX
f16s <- taxa16S
rownames(f16s) <- f16s[,1]
f16s[,2] <- apply(f16s[,1:2], 1, paste, collapse=' ')
f16s <- gsub(' ', '_', f16s)
df <- data.frame(ASV=apply(f16s[rownames(t16s),-1], 1, paste, collapse=';'), t16s)
write.table(df,
            sprintf("%s/FAPROTAX_Input.tsv", dir$tabledir), 
            row.names=FALSE, quote=FALSE, sep="\t")


## -- For FunGUILD
tits <- t(dflist[[2]])

taxaITSsub <- taxaITS[,-1]
initial <- c("k__",'p__','c__','o__','f__','g__','s__')
for(i in 1:ncol(taxaITSsub)){
    taxaITSsub[,i] <- paste(initial[i], taxaITSsub[,i],sep="")
}
taxainfo <- apply(taxaITSsub[rownames(tits), ], 1, paste, collapse=";")
funguild <- cbind(ID=rownames(tits), tits, taxonomy=taxainfo)

write.table(funguild,
            sprintf("%s/funguildInput.txt", dir$tabledir), 
            row.names=FALSE, quote=FALSE, sep="\t")

