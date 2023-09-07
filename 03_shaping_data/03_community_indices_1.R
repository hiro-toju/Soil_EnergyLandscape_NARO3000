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
load.lib(c('ggplot2', 'tidyr', 'vegan'))

# -- Create directory to save
dir <- make.dir('03_shaping_data/03_output')

## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# -- Load data table

## -- Converted
abs16s <- as.data.frame(readRDS('Table/readConvert_Prok.rds'))
absits <- as.data.frame(readRDS('Table/readConvert_Fng.rds')) 
abslist <- list(abs16s, absits)

rarefied16s <- readRDS('Table/rrarefy_16S.rds')
rarefiedits <- readRDS('Table/rrarefy_ITS.rds')
rarefied <- list(rarefied16s, rarefiedits)

saveRDS(rarefied, 'Table/rrarefied_Matrix.rds')
sml <- readRDS("Table/0301_sml.rds") ; dim(sml)

################################################################

relative <- lapply(rarefied, function(x){ x/rowSums(x) })
lapply(relative, dim)
## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
## -- Alpha diversity
alpha <- lapply(rarefied, function(x){
				as.data.frame(cbind(richness = apply(x, 1, function(x){sum(x > 0)}),
							  		simpson  = diversity(x, index ='simpson'),
			 				  		shannon  = diversity(x, index ='shannon')
							  )	)	})

index <- cbind(bac=alpha[[1]][rownames(sml),], 
				   bac.total.abundance=rowSums(abs16s[rownames(sml),]),
				   fng=alpha[[2]][rownames(sml),], 
				   fng.total.abundance=rowSums(absits[rownames(sml),]) )
				   
ASV.ID <- rownames(sml)
				
indices <- cbind(index, BbyFabundance=log(index[,c('bac.total.abundance')])-log(index[,c('fng.total.abundance')]),
						BbyFrichness=(index[,c('bac.richness')])/(index[,c('fng.richness')]))

rownames(indices) <- ASV.ID

pdf(sprintf('%s/bacteria_fungi_rrarefy.pdf', dir$figdir), w=15, h=15)
pairs(indices); dev.off()

gg <- ggplot(indices, aes(x= bac.shannon, y= fng.shannon))+
	  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
	  geom_density_2d(color='grey90', alpha=0.6)+
	  geom_point(alpha=0.5, color='orange', size=1)+
	  scale_fill_continuous(type = "viridis")+
	  scale_x_continuous(expand = c(0, 0)) +
	  scale_y_continuous(expand = c(0, 0)) +
	  theme(legend.position='none',
	  		axis.title=element_text(size=15) )+
	  labs(x='Bacteria alpha diveristy', y='Fungi alpha diveristy')
	  
	  
ggsave(plot=gg, filename=sprintf('%s/bacteria_fungi_alphaDiv_rrarefy.pdf', dir$figdir), w=7, h=7)

saveRDS(indices, 'Table/0303_communityIndices_rrarefy.rds')
saveRDS(relative, 'Table/0303_rarefiedMatrix_rrarefy.rds')
#saveRDS(relative, 'Table/0303_rarefiedMatrix.rds')



write.table(cbind(ASV.ID, indices), file=sprintf('%s/0303_communityIndices_rrarefy.txt', dir$tabledir), quote=F, sep="\t", row.names=F)

############################################################################