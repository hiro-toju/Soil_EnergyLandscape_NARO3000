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

l=${l}

setwd('/user4/kyoto1/hiroakif/aptmp/NARO')
############################################################################

# -- Create directory to save
outdir="${outdir}"
dir <- make.dir( sprintf('%s/02_output', outdir))

# -- Load data table
result01 <- list.files(  sprintf("%s/01_output/RDS/", outdir), full.names = TRUE)

sml <- readRDS("Table/0301_sml.rds") 

############################################################################

tL <- c('Family', 'Order')
plun=0.1
scale2 <- function(X){ (X-min(X,na.rm=TRUE))/(max(X,na.rm=TRUE)-min(X,na.rm=TRUE)) }
thread=${ncpu}/2

ConvDec <- function(x){sep <- round(length(x)/20)
line <- c()
for(i in 1:sep){
    y <- x[-length(x)]
    bi <- na.omit(y[((1+20*(i-1)):(20*i))])
    line <- c(line, strtoi( paste0(bi, collapse=''), base=2))
    
}

ssid=paste(line, collapse='-')
c(ssid) }

############################################################################

for(tl in rev(tL)){ #tl='Family'
    
    res <- readRDS(result01[grep(tl, result01)])
    
    fit <- res[[1]]
    ocdata <- res[[2]]
    
    envstd <- res[[3]]    
    
    stab <- ssmat <- c()
    for(i in rownames(ocdata)){ #
        #i=rownames(ocdata)[l]
        
        ## |||||||||||||||||||||||||||||||||||||| ##
        com <- as.numeric(ocdata[i,])
        ee <- envstd[i,]
        
        ## |||||||||||||||||||||||||||||||||||||| ##
        
        alpha <- fit[,grep('^g\\.', colnames(fit))]%*%as.matrix(ee) + fit[,'h']
        beta <- fit[,grep('^J\\.', colnames(fit))]
        
        ## |||||||||||||||||||||||||||||||||||||| ##
        minsets = rELA:::SSestimate(alpha, beta, itr = 20000)
        minsets = unique(minsets)
        uniqueSS <- data.frame(unique(minsets[, -ncol(minsets)]))
        
        ## |||||||||||||||||||||||||||||||||||||| ##        
        sampleSS <- rELA:::SteepestDescent_cpp(com, alpha=alpha, beta=beta)
        ssid <- ConvDec(sampleSS)
        
        ssenergy <- sampleSS[, ncol(sampleSS)]
        energy <- rELA:::cEnergy(com, alpha = alpha, beta = beta)
        energy.gap <- energy - ssenergy
        
        res <- data.frame(StableState.id = ssid, SSenergy = ssenergy, 
                          energy=energy, 
                          energy.gap = energy.gap, nSS=nrow(uniqueSS))    				
        stab <- rbind(stab, res)
        ssmat <- rbind(ssmat, sampleSS)
    }
    
    rownames(stab) <- rownames(ocdata)
    saveRDS(stab, sprintf('%s/0602_ELAstability_%s.rds', dir$rdsdir, tl))
    
    rownames(ssmat) <- rownames(ocdata)
    colnames(ssmat) <- c(colnames(ocdata), 'energy')
    saveRDS(ssmat, sprintf('%s/0602_SScomm_%s.rds', dir$rdsdir, tl))
    
}
## ====================================================== ##
