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

## -- For super computer system

# -- Load library and functions
#　lapply(lib, function(x){ install.packages(lib, lib=lib.at)})
lib <- c('doParallel', 'rELA', 'AnalysisHelper')
lib.at <- '/user4/kyoto1/hiroakif/R/x86_64-pc-linux-gnu-library/4.1'

invisible(lapply(lib, function(x){ library(package = x, character.only=TRUE, lib=lib.at)}) )
invisible(sapply(lib, function(x) cat(sprintf('%s %s\n',x, packageVersion(x, lib=lib.at)) ) ))

setwd('/user4/kyoto1/hiroakif/aptmp/NARO')
############################################################################

# -- Create directory to save
outdir="${outdir}"
dir <- make.dir( sprintf('%s/01_output', outdir))

# -- Load data table
#read16s <- as.data.frame(readRDS('Table/read_Prok.rds'))
readdf <- as.data.frame(readRDS('Table/rrarefy_ITS.rds')) 

#commonrow <- intersect(rownames(read16s), rownames(readdf))

#dflist <- list(read16s[commonrow, colSums(read16s)>0], readdf[commonrow, colSums(readdf)>0])

sml <- readRDS("Table/0301_sml.rds") 

taxa <- readRDS('Table/taxonomyList_Fng.rds')

############################################################################

targetENV <- c('pH_dry_soil', 'CNratio', 'available_P','EC_electric_conductivity')

tL <- c('Family', 'Order')

min.filtth=0.05
max.filtth=0.9
itr=100000

thread=8


############################################################################

# -- Filtering OTUs for analysis of community.
filt <- function(x, th=10, nsample=nrow(x),
                 minth=0.05, maxth=1){
    freq <- colSums(x>th)/nsample
    cat( sprintf('ASVs frequency range : minimum is %s, maximum is %s\n', 
                 round(min(freq), 2), round(max(freq), 2)))
    cat(sprintf('Discard ASVs in which appeared less than %s samples\n', round(nsample*minth)))
    y <- x[, freq>minth & freq<maxth]
    z <- y[rowSums(y)>0, ]
    cat( sprintf('Remained sample/asvs is %s/%s\n', nrow(z), ncol(z)))
    return(z) }

scale2 <- function(X){ (X-min(X,na.rm=TRUE))/(max(X,na.rm=TRUE)-min(X,na.rm=TRUE)) }

################################################################

## -- Environement samples
smltar <-(sml[,targetENV])
envrow <- rownames(na.omit(smltar))
#envrow <- rownames(sml)[!is.na((sml[,targetENV]))]

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
## -- Non-NA samples
commonrow <- intersect(rownames(readdf), envrow)

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##

for(tl in tL){ #tl='Genus'
    
    ## -- Merge matrix and filtering and binarize
    dfcommon <- Taxa.mat(readdf, taxa, tl)
    dfcommon <- dfcommon[commonrow,]
    
    if(any(colnames(dfcommon)=='Unidentified')){
        dfcommon <- dfcommon[,which(colnames(dfcommon)!='Unidentified')]
    }
    
    oclist <- filt(dfcommon, minth=min.filtth, th=0, maxth=max.filtth)
    commonrow <- intersect(rownames(oclist), envrow)
    
    ocdata <- oclist[commonrow,]
    ocdata[ocdata>0] <- 1
    
    ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
    ## -- Standarization
    env <- sml[, targetENV]
    envstd <- apply(as.matrix(env), 2, scale2)#[commonrow,]
    rownames(envstd) <- rownames(sml)
    
    envstd <- envstd[commonrow,]
    
    ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
    ## Parameter fitting (This proceess takes a long time)
    fittingRes <- runSAparallel(as.matrix(ocdata), env=envstd,
                                rep = 20, max.itr = itr, thread=thread)
    
    saveRDS(list(fittingRes=fittingRes, community=ocdata, env=envstd), sprintf('%s/0601_ELAfitting_%s.rds', dir$rdsdir, tl) )
    
}


## ====================================================== ##
