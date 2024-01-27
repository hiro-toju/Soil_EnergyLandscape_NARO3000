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

## -- Create directory to save
dir <- make.dir('05_Community_structure/M2_02_output')

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

col1 <- brewer.pal(7, "Set1")
col2 <- brewer.pal(6, "Set2")
col3 <- brewer.pal(11, "Set3")
col4 <- brewer.pal(8, "Pastel1")
col5 <- brewer.pal(3, "Pastel1")

colvec <- c(col1, col2, col3, col4, col5)

############################################################################

## -- Merge level
taxalevel <- 'Genus'
n.core <- parallel::detectCores()

############################################################################

## -- Barplot
mergelist <- lapply(dflist, Taxa.mat, taxa=taxa, taxaLabel=taxalevel)

cols <- lapply(mergelist, makeColPalette, color= colvec,
               specificName = c('Bac_Unidentified', 'Arc_Unidentified','Fng_Unidentified'), 
               specificColor = c('grey90'))

rellist <- c()
for(i in 1:2){
    
    dftmp <- mergelist[[i]]/rowSums(mergelist[[i]])
    coltmp <- cols[[i]]
    rellist[[i]] <- dftmp
    if(any(coltmp=='grey30')){
        rellist[[i]] <- cbind(dftmp[ ,names(coltmp[-which(coltmp=='grey30')]) ], 
                              'Others'=rowSums(dftmp[ ,names(coltmp[which(coltmp=='grey30')]) ]))
    }
    
}

colsub <- lapply(cols, function(x){ #x=cols[[1]]
    
    if(any(x=='grey30')) { res <- c(x[-which(x=='grey30')], 'Others'='grey30')
    }else{res <- x}
    return(res)
})

lflist <- lapply(rellist, function(x){ #x=rellist[[1]]
    
    tmp <- cbind(sml, x[rownames(x),])
    lf <- gather(tmp, key, value, -c(1:ncol(sml)))
     
    return(lf)
})

ggbase <- function(data=data, ..., base_size=10){
            ggplot(data, ...)+
            scale_x_discrete(expand=c(0,0))+
            scale_y_continuous(expand=c(0,0))+
            theme_minimal(base_size=base_size)+
            theme(panel.border=element_rect(fill=NA, size=0.6)) }

############################################################################
# without legends

gglist <- lapply(1:2, function(x){ #x=1
        
	data=lflist[[x]]
	dataagg <- aggregate(data$value, by=list(env=data[,'accession_number'], key=data$key), sum, na.rm=TRUE)

    pal=colsub[[x]]
    dataagg $key <- factor(dataagg $key, levels=rev(names(pal)))
    dataagg $env <- factor(dataagg $env, levels=sort(unique(dataagg $env)))
    
ggg <- ggbase(dataagg, aes(x= env, y=x))+
            geom_bar( aes(fill=key), 
                      color=NA, width=1,size=0.1, 
                      stat='identity', position='fill')+
            scale_fill_manual(values=pal)+
                      labs(x='', y='Proportion of sequencing reads', fill=taxalevel)+
            guides(fill=guide_legend(title='',ncol=1, reverse=T)) + theme(legend.position = "none")

    return(ggg)
})
    
ggsave(plot=plot_grid(plotlist=gglist, ncol=1, align='v'),
       filename=sprintf('%s/barplot_Genus_rrarefy.pdf', dir$figdir), w=6, h=5)

############################################################################
# with legends

gglist <- lapply(1:2, function(x){ #x=1
        
	data=lflist[[x]]
	dataagg <- aggregate(data$value, by=list(env=data[,'accession_number'], key=data$key), sum, na.rm=TRUE)

    pal=colsub[[x]]
    dataagg $key <- factor(dataagg $key, levels=rev(names(pal)))
    dataagg $env <- factor(dataagg $env, levels=sort(unique(dataagg $env)))
    
ggg <- ggbase(dataagg, aes(x= env, y=x))+
            geom_bar( aes(fill=key), 
                      color=NA, width=1,size=0.1, 
                      stat='identity', position='fill')+
            scale_fill_manual(values=pal)+
                      labs(x='', y='Proportion of sequencing reads', fill=taxalevel)+
            guides(fill=guide_legend(title='',ncol=1, reverse=T)) + theme(legend.position = "right")

    return(ggg)
})
    
ggsave(plot=gglist[[1]],
       filename=sprintf('%s/barplot_Genus_rrarefy_pro_caption.pdf', dir$figdir), w=6, h=10)
       
ggsave(plot=gglist[[2]],
       filename=sprintf('%s/barplot_Genus_rrarefy_fng_caption.pdf', dir$figdir), w=6, h=10)


